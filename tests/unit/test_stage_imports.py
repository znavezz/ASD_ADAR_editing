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
