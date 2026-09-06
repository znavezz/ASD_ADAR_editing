"""The documented entry points must actually start.

`scripts/run_pipeline.py` was unrunnable for the whole of the directory
reorganisation: it imports `helpers`, which moved into `pipeline/`, and nothing
put that directory on `sys.path` outside pytest's own conftest. Parse checks and
the rest of the suite all passed, because nothing invoked the script. These tests
invoke it.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _run(*args):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_pipeline.py"), *args],
        capture_output=True, text=True, cwd=ROOT, timeout=180,
    )


def test_runner_starts_and_reports_its_options():
    """--help must succeed, which requires every module-level import to resolve."""
    r = _run("--help")
    assert r.returncode == 0, f"entry point does not start:\n{r.stderr}"
    assert "--env-file" in r.stdout
    assert "--start_from" in r.stdout


def test_runner_refuses_a_missing_env_file():
    """The env file names the database written to, so a bad one must stop the run."""
    r = _run("--env-file", "does-not-exist.env")
    assert r.returncode != 0
    assert "env file not found" in (r.stderr + r.stdout)


def test_reset_db_script_is_valid_shell():
    r = subprocess.run(["bash", "-n", str(ROOT / "reset_db.sh")],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr


def test_reset_db_exports_its_environment():
    """The env file must reach child processes, not just the shell.

    `docker compose` does not see unexported shell variables; it reads ./.env
    itself. Sourcing .env.test without exporting truncated the test database while
    resolving COMPOSE_PROJECT_NAME to the production project, so `docker compose
    down` targeted the wrong stack. `set -a` around the source is what prevents it.
    """
    text = (ROOT / "reset_db.sh").read_text()
    source_at = text.index('source "$ENV_FILE"')
    assert "set -a" in text[:source_at], "env file is sourced without being exported"
    assert "set +a" in text[source_at:], "allexport is never turned back off"
