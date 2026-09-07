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


def test_shell_entry_points_are_executable():
    """A script the README invokes as ./name must carry the executable bit in git.

    reset_db.sh was committed 100644, so the documented `./reset_db.sh` failed with
    Permission denied in every fresh clone. The bit lives in the index, not just on
    disk, so chmod alone does not fix it for anyone else.
    """
    import stat
    for rel in ("reset_db.sh", "figures/Lit_Search/06_mode_b_audit_runner_ver2.sh"):
        mode = (ROOT / rel).stat().st_mode
        assert mode & stat.S_IXUSR, f"{rel} is not executable"


def test_reset_refuses_the_published_database():
    """Naming the published database must not be enough to delete it.

    --before 1 is a full teardown that removes the PostgreSQL data directory. The
    default target is the scratch stack, so this asks for the real one by name and
    asserts that asking is still refused without the opt-in.
    """
    r = subprocess.run(
        ["bash", str(ROOT / "reset_db.sh"), "--before", "1", "--env-file", ".env"],
        capture_output=True, text=True, cwd=ROOT, timeout=120)
    assert r.returncode != 0, "a full reset of the published database was not refused"
    assert "published database" in r.stderr


def test_destructive_entry_points_default_to_the_scratch_stack():
    """Both the reset and the pipeline write. Their default target must be the test
    stack, not the database that produced the manuscript."""
    assert 'ASD_ENV_FILE:-.env.test' in (ROOT / "reset_db.sh").read_text()
    assert '"ASD_ENV_FILE", ".env.test"' in (ROOT / "scripts" / "run_pipeline.py").read_text()


def test_pipeline_refuses_to_write_to_the_published_database():
    r = _run("--env-file", ".env")
    assert r.returncode != 0
    assert "published database" in (r.stderr + r.stdout)


def test_reset_wait_loop_survives_set_e():
    """`((i++))` returns the pre-increment value, so at i=0 it exits 1 and `set -e`
    kills the script. The container-start path died on its first iteration and had
    therefore never worked."""
    code = [l for l in (ROOT / "reset_db.sh").read_text().splitlines()
            if not l.lstrip().startswith("#")]
    offenders = [l for l in code if "((i++))" in l or "((j++))" in l]
    assert not offenders, f"post-increment in a set -e script aborts at 0: {offenders}"
