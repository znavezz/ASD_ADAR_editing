#!/usr/bin/env python
"""Capture and verify the published-number oracle.

The paper's headline numbers were produced by hand-running the GraphQL operations in
``queries.txt`` against Hasura. That file is inert text, so it can drift from the numbers
recorded in its own comments -- and it had: the off-target filter in
``NonG2A_StopGained_Pass_NoDeleteriousBystander`` was commented out, making the query return
404 instead of the 315 its comment claimed.

This tool makes the oracle executable:

* ``capture``  runs every named operation and writes ``tests/golden/oracle/counts.json``
  plus ``tests/golden/oracle/artifacts.sha256`` (checksums of the gitignored Output/ CSVs).
* ``verify``   re-runs everything and fails if any count moved, or if any ``# expect:``
  annotation in queries.txt disagrees with the live result.

Annotate an expected value inside a query block with either form::

    # expect: 2421                 (when the query returns exactly one count)
    # expect TotalVariants: 329279 (when it returns several aliases)

Usage::

    python scripts/verify_counts.py capture
    python scripts/verify_counts.py verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
QUERIES_FILE = REPO_ROOT / "queries" / "queries.txt"
ORACLE_DIR = REPO_ROOT / "tests" / "golden" / "oracle"
COUNTS_FILE = ORACLE_DIR / "counts.json"
ARTIFACTS_FILE = ORACLE_DIR / "artifacts.sha256"

#: Output/ artifacts that back published figures and supplementary tables. These are
#: gitignored (they are large and regenerable), so we commit their checksums instead.
#:
#: The figure source-data files live under ``Output/Figures_v1_2026-05-20/``, whose
#: directory names carry the figure numbers *as of that run* -- Figure2 there is the
#: CADD/rescue-SIFT scatter, which the 2026-08-16 split renumbered to Figure 3. The
#: directory is deliberately not renumbered: it is the frozen record of the tree that
#: produced the published values, and renaming it to today's numbering would make a
#: checksum claim about a file that run never wrote. Current-generation output goes to
#: ``Output/Figures/Figure_<n>/`` instead, so the two never collide.
GOLDEN_ARTIFACTS = [
    "Output/Intermediate_tables/00_db_snapshot/df_cds.csv",
    # Hashed over its content, not its bytes -- see hash_file().
    "Output/Intermediate_tables/00_db_snapshot/manifest.json",
    "Output/Figures_v1_2026-05-20/Figure1/Source_data_Figure_1A_funnel.csv",
    "Output/Figures_v1_2026-05-20/Figure1/Source_data_Figure_1B_mismatch.csv",
    "Output/Figures_v1_2026-05-20/Figure1/Source_data_Figure_1C_clinvar_donut.csv",
    "Output/Figures_v1_2026-05-20/Figure1/Source_data_Figure_1D_offtargets.csv",
    "Output/Figures_v1_2026-05-20/Figure1/Totals_Figure_1B_mismatch.csv",
    "Output/Figures_v1_2026-05-20/Figure2/Source_data_Figure_2B_cadd_sift.csv",
    "Output/Figures_v1_2026-05-20/Figure2/Supp_Table_F2B_Figure_2B_cadd_sift_quadrants.csv",
    "Output/Figures_v1_2026-05-20/Figure3/Source_data_Figure_3_brain_heatmap.csv",
    "Output/Figures_v1_2026-05-20/Figure3/Source_data_Figure_3_treatability_bar.csv",
    "Output/Figures_v1_2026-05-20/Figure3/Supp_Table_F3B_Figure_3_treatability_bar.csv",
    "Output/TableS1_Fixable_variants.csv",
    "Output/TableS2_Rescue_variants.csv",
    "Output/TableS3_Missense_optimization_variants.csv",
    "Output/Intermediate_tables/panels/MissenseOptimization_full_flat.csv",
]


# --------------------------------------------------------------------------------------
# Parsing queries.txt
# --------------------------------------------------------------------------------------

_QUERY_START = re.compile(r"^\s*(?:query|mutation)\s+(\w+)\s*[({]")
_EXPECT = re.compile(r"#\s*expect\s*(?:\s(\w+))?\s*:\s*([\d,_]+)\s*$")


@dataclass
class Query:
    """One named GraphQL operation lifted out of queries.txt."""

    name: str
    text: str
    line: int
    expect: dict[str | None, int] = field(default_factory=dict)
    has_variables: bool = False

    def check(self, counts: dict[str, int]) -> list[str]:
        """Return a list of human-readable problems; empty means the query agrees with
        its own ``# expect:`` annotations."""
        problems: list[str] = []
        for alias, expected in self.expect.items():
            if alias is None:
                if len(counts) != 1:
                    problems.append(
                        f"{self.name} (queries.txt:{self.line}): bare '# expect: {expected}' is "
                        f"ambiguous -- the query returned {len(counts)} counts "
                        f"({', '.join(sorted(counts)) or 'none'}). Use '# expect <alias>: N'."
                    )
                    continue
                actual = next(iter(counts.values()))
                key = next(iter(counts))
            else:
                key = _resolve_alias(alias, counts)
                if key is None:
                    problems.append(
                        f"{self.name} (queries.txt:{self.line}): '# expect {alias}: {expected}' "
                        f"names an alias that is not in the result "
                        f"({', '.join(sorted(counts)) or 'none'})."
                    )
                    continue
                actual = counts[key]
            if actual != expected:
                problems.append(
                    f"{self.name}.{key} (queries.txt:{self.line}): expected {expected:,}, "
                    f"got {actual:,}."
                )
        return problems


def _resolve_alias(alias: str, counts: dict[str, int]) -> str | None:
    if alias in counts:
        return alias
    for key in counts:
        if key.split(".")[0] == alias:
            return key
    return None


def _strip_noise(line: str) -> str:
    """Blank out string literals and the comment tail so brace counting is accurate.

    Both matter in practice: queries.txt filters contain braces inside `_eq: "..."`
    values, and a commented-out filter line carries unbalanced braces.
    """
    out: list[str] = []
    in_string = False
    escaped = False
    for ch in line:
        if in_string:
            out.append(" ")
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            continue
        if ch == '"':
            in_string = True
            out.append(" ")
            continue
        if ch == "#":
            break
        out.append(ch)
    return "".join(out)


def parse_queries(text: str) -> list[Query]:
    """Split a GraphQL document into its named operations.

    Brace depth is tracked with string literals and comments blanked out, so neither a
    ``"a{b"`` value nor a commented-out filter line can end a block early.
    """
    lines = text.split("\n")
    queries: list[Query] = []
    i = 0
    while i < len(lines):
        match = _QUERY_START.match(lines[i])
        if not match:
            i += 1
            continue

        start = i
        depth = 0
        opened = False
        block: list[str] = []
        while i < len(lines):
            raw = lines[i]
            block.append(raw)
            clean = _strip_noise(raw)
            depth += clean.count("{") - clean.count("}")
            if "{" in clean:
                opened = True
            i += 1
            if opened and depth <= 0:
                break

        body = "\n".join(block)
        expect: dict[str | None, int] = {}
        for raw in block:
            found = _EXPECT.search(raw)
            if found:
                alias, value = found.groups()
                expect[alias] = int(value.replace(",", "").replace("_", ""))

        queries.append(
            Query(
                name=match.group(1),
                text=body,
                line=start + 1,
                expect=expect,
                has_variables="(" in lines[start].split(match.group(1), 1)[1].split("{")[0],
            )
        )
    return queries


# --------------------------------------------------------------------------------------
# Result reduction
# --------------------------------------------------------------------------------------


def extract_counts(data: object, prefix: str = "") -> dict[str, int]:
    """Reduce a GraphQL response to ``{path: count}``.

    An ``{"aggregate": {"count": N}}`` node contributes its parent path; a list
    contributes its length. Everything else is walked through.
    """
    counts: dict[str, int] = {}
    if isinstance(data, dict):
        agg = data.get("aggregate")
        if isinstance(agg, dict) and isinstance(agg.get("count"), int):
            counts[prefix] = agg["count"]
            return counts
        for key, value in data.items():
            path = f"{prefix}.{key}" if prefix else key
            counts.update(extract_counts(value, path))
    elif isinstance(data, list):
        if prefix:
            counts[prefix] = len(data)
    return counts


# --------------------------------------------------------------------------------------
# Hasura
# --------------------------------------------------------------------------------------


def load_env(env_file: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not env_file.exists():
        return env
    for line in env_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        env[key.strip()] = value.strip().strip('"').strip("'")
    return env


class Hasura:
    def __init__(self, url: str, secret: str, timeout: int = 300):
        self.url = url
        self.secret = secret
        self.timeout = timeout

    def run(self, query: str, variables: dict | None = None) -> dict:
        import requests

        response = requests.post(
            self.url,
            headers={"x-hasura-admin-secret": self.secret, "Content-Type": "application/json"},
            json={"query": query, "variables": variables or {}},
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()
        if "errors" in payload:
            raise RuntimeError(f"GraphQL error: {json.dumps(payload['errors'])[:500]}")
        return payload.get("data") or {}


def connect(env_file: Path) -> Hasura:
    env = load_env(env_file)
    port = os.environ.get("GRAPHQL_PORT") or env.get("GRAPHQL_PORT", "8789")
    host = os.environ.get("SERVER_HOST") or env.get("SERVER_HOST", "127.0.0.1")
    secret = os.environ.get("HASURA_ADMIN_SECRET") or env.get("HASURA_ADMIN_SECRET", "")
    if not secret:
        raise SystemExit(f"No HASURA_ADMIN_SECRET found in {env_file} or the environment.")
    return Hasura(f"http://{host}:{port}/v1/graphql", secret)


# --------------------------------------------------------------------------------------
# Artifacts
# --------------------------------------------------------------------------------------


#: Manifest keys that describe *this run* rather than the snapshot: a wall-clock
#: timestamp and the endpoint it happened to be pulled from. Hashing them would
#: make every honest re-run look like a regression; dropping the manifest from the
#: pinned set entirely would stop the facts it *does* carry -- row counts, the
#: consequence breakdown, the query's own sha1 -- from being checked at all. So the
#: checksum covers the content and skips the circumstances.
_MANIFEST_RUN_CONTEXT = ("pulled_at", "hasura_url")


def hash_file(path: Path) -> tuple[str, int]:
    """Return (sha256, line_count) for a file, streamed.

    ``manifest.json`` is hashed over a canonical projection of its content rather
    than its bytes, so the digest states what the snapshot contains instead of
    when it was taken.
    """
    if path.name == "manifest.json":
        raw = path.read_text()
        content = {
            k: v for k, v in json.loads(raw).items()
            if k not in _MANIFEST_RUN_CONTEXT
        }
        payload = json.dumps(content, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(payload).hexdigest(), raw.count("\n")

    digest = hashlib.sha256()
    lines = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
            lines += chunk.count(b"\n")
    return digest.hexdigest(), lines


def capture_artifacts() -> list[str]:
    rows = [
        "# sha256  lines  path",
        "# Checksums of the Output/ artifacts that back the published figures and tables.",
        "# Output/ is gitignored, so these checksums are the committed record.",
        "# Regenerate with: python scripts/verify_counts.py capture",
    ]
    for rel in GOLDEN_ARTIFACTS:
        path = REPO_ROOT / rel
        if not path.exists():
            rows.append(f"MISSING  -  {rel}")
            continue
        digest, lines = hash_file(path)
        rows.append(f"{digest}  {lines}  {rel}")
    return rows


# --------------------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------------------


def run_all(hasura: Hasura, queries: list[Query]) -> tuple[dict, list[str]]:
    results: dict[str, dict] = {}
    problems: list[str] = []
    for query in queries:
        if query.has_variables:
            results[query.name] = {
                "skipped": "requires variables",
                "line": query.line,
                "query_sha1": hashlib.sha1(query.text.encode()).hexdigest(),
            }
            print(f"  {query.name:<48} skipped (requires variables)")
            continue
        try:
            data = hasura.run(query.text)
        except Exception as exc:  # noqa: BLE001 - report and continue over all queries
            results[query.name] = {"error": str(exc)[:300], "line": query.line}
            problems.append(f"{query.name} (queries.txt:{query.line}): {exc}")
            print(f"  {query.name:<48} ERROR")
            continue

        counts = extract_counts(data)
        results[query.name] = {
            "line": query.line,
            "query_sha1": hashlib.sha1(query.text.encode()).hexdigest(),
            "counts": counts,
            "expect": {str(k): v for k, v in query.expect.items()},
        }
        problems.extend(query.check(counts))
        rendered = ", ".join(f"{k}={v:,}" for k, v in counts.items())
        print(f"  {query.name:<48} {rendered}")
    return results, problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("mode", choices=["capture", "verify"])
    parser.add_argument("--env-file", type=Path, default=REPO_ROOT / ".env")
    parser.add_argument(
        "--skip-artifacts", action="store_true", help="Only handle query counts."
    )
    args = parser.parse_args(argv)

    queries = parse_queries(QUERIES_FILE.read_text())
    print(f"Parsed {len(queries)} operations from {QUERIES_FILE.name}\n")

    hasura = connect(args.env_file)
    results, problems = run_all(hasura, queries)

    payload = {
        "hasura_url": hasura.url,
        "queries_sha256": hashlib.sha256(QUERIES_FILE.read_bytes()).hexdigest(),
        "results": results,
    }

    if args.mode == "capture":
        ORACLE_DIR.mkdir(parents=True, exist_ok=True)
        COUNTS_FILE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"\nWrote {COUNTS_FILE.relative_to(REPO_ROOT)}")
        if not args.skip_artifacts:
            ARTIFACTS_FILE.write_text("\n".join(capture_artifacts()) + "\n")
            print(f"Wrote {ARTIFACTS_FILE.relative_to(REPO_ROOT)}")
    else:
        problems.extend(_verify_against_baseline(results))
        if not args.skip_artifacts:
            problems.extend(_verify_artifacts())

    if problems:
        print(f"\n{len(problems)} problem(s):")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print("\nAll counts and annotations agree.")
    return 0


def _verify_against_baseline(results: dict) -> list[str]:
    if not COUNTS_FILE.exists():
        return [f"No baseline at {COUNTS_FILE.relative_to(REPO_ROOT)}; run `capture` first."]
    baseline = json.loads(COUNTS_FILE.read_text())["results"]
    problems = []
    for name, current in results.items():
        want = baseline.get(name, {}).get("counts")
        got = current.get("counts")
        if want is None or got is None:
            continue
        for key, expected in want.items():
            actual = got.get(key)
            if actual != expected:
                problems.append(
                    f"{name}.{key}: baseline {expected:,}, now "
                    f"{'missing' if actual is None else format(actual, ',')}."
                )
    return problems


def _verify_artifacts() -> list[str]:
    if not ARTIFACTS_FILE.exists():
        return []
    problems = []
    for line in ARTIFACTS_FILE.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        digest, lines, rel = line.split(None, 2)
        path = REPO_ROOT / rel
        if digest == "MISSING":
            continue
        if not path.exists():
            problems.append(f"{rel}: recorded but now missing.")
            continue
        actual_digest, actual_lines = hash_file(path)
        if actual_digest != digest:
            problems.append(
                f"{rel}: sha256 changed ({lines} -> {actual_lines} lines)."
            )
    return problems


if __name__ == "__main__":
    sys.exit(main())
