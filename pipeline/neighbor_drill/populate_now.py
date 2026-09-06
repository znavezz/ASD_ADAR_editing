#!/usr/bin/env python3
"""
THROWAWAY backfill — populate `neighbor_edits` + Improve guides/bystanders on the ALREADY-populated
DB, without a full pipeline re-run.

The canonical logic lives in phase3 (via `neighbor_drill.*_pipeline`); this just calls the same
pipeline functions against the current DB. Delete once phase3's embedding is verified end-to-end.

    python neighbor_drill/populate_now.py --env-file .env.test [--limit N] [--dry-run]
                                          [--only neighbor_edits|improve_guides]
"""

import sys
import argparse
import logging
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
from dotenv import load_dotenv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

ap = argparse.ArgumentParser()
ap.add_argument("--env-file", default=".env")
ap.add_argument("--limit", type=int, default=None)
ap.add_argument("--dry-run", action="store_true")
ap.add_argument("--only", choices=["neighbor_edits", "improve_guides"], default=None)
args = ap.parse_args()

env_path = Path(args.env_file)
if not env_path.is_absolute():
    env_path = _ROOT / env_path
if not env_path.exists():
    raise RuntimeError(f"env file not found: {env_path}")
load_dotenv(env_path, override=True)
logging.info(f"Loaded env {env_path.name}")

from neighbor_drill.neighbor_edits_pipeline import populate_neighbor_edits
from neighbor_drill.improve_guides_pipeline import populate_improve_guides

if args.only in (None, "neighbor_edits"):
    logging.info("=== populate neighbor_edits ===")
    populate_neighbor_edits(limit=args.limit, dry_run=args.dry_run)
if args.only in (None, "improve_guides"):
    logging.info("=== populate Improve guides + bystanders ===")
    populate_improve_guides(limit=args.limit, dry_run=args.dry_run)
