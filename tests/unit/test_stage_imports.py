"""The stages and the runner must bind shared names the same way.

The stages used to inherit names from the runner's scope via `exec()`. When they
were given explicit imports, `phase2_db.py` acquired `import datetime` where the
runner had `from datetime import datetime`, so every `datetime.utcnow()` call
became an AttributeError on the module. Nothing caught it: the file parsed, the
suite passed, and the failure surfaced only two hours into a full pipeline run.

A name bound as a module in one file and as a class in another is the signature
of that mistake, so this test compares the bindings directly.
"""
import ast
import collections
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = [ROOT / "scripts" / "run_pipeline.py"] + sorted((ROOT / "pipeline").glob("phase*.py"))


def _bindings():
    """name -> {file: binding-shape} across the runner and every stage."""
    out = collections.defaultdict(dict)
    for f in FILES:
        for node in ast.walk(ast.parse(f.read_text())):
            if isinstance(node, ast.Import):
                for a in node.names:
                    out[a.asname or a.name.split(".")[0]][f.name] = f"module:{a.name}"
            elif isinstance(node, ast.ImportFrom):
                for a in node.names:
                    out[a.asname or a.name][f.name] = f"from {node.module} import {a.name}"
    return out


def test_shared_names_are_bound_consistently():
    conflicts = {
        name: per for name, per in _bindings().items() if len(set(per.values())) > 1
    }
    assert not conflicts, (
        "these names are bound differently across the runner and the stages, so a "
        f"call valid in one file raises in another: {conflicts}"
    )


def test_phase2_binds_datetime_as_the_class():
    """The specific regression: phase2_db.py calls datetime.utcnow()."""
    src = (ROOT / "pipeline" / "phase2_db.py").read_text()
    assert "from datetime import datetime" in src
    assert "\nimport datetime\n" not in src


def test_data_paths_resolve_against_the_repository_root():
    """neighbor_drill resolves relative data paths like "Resources/...".

    `_ROOT` there is pipeline/ - correct for the sys.path insert, since helpers.py
    lives beside it - but data paths are relative to the repository, one level
    further up. Using one variable for both resolved Resources/... to
    pipeline/Resources/... once neighbor_drill moved a level deeper, and the
    pipeline died on it nine hours into a run.
    """
    import sys
    sys.path.insert(0, str(ROOT / "pipeline"))
    from neighbor_drill.neighbor_edits_pipeline import _resolve as resolve_edits
    from neighbor_drill.improve_guides_pipeline import _resolve as resolve_guides

    for resolve in (resolve_edits, resolve_guides):
        got = resolve("Resources/Homo_sapiens.GRCh37.87.gtf.gz")
        assert got == ROOT / "Resources" / "Homo_sapiens.GRCh37.87.gtf.gz", got
        assert "pipeline/Resources" not in str(got)
