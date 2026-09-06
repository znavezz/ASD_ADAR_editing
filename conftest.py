"""Root conftest: put the repository root and pipeline/ on sys.path so tests can
import `pipeline`, `scripts` and `neighbor_drill` without an install step.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for entry in (ROOT, ROOT / "pipeline"):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))
