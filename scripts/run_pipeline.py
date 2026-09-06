## PreProcess
import sys
from pathlib import Path as _Path

# The stages live in pipeline/ and import each other by bare name (`helpers`,
# `neighbor_drill`). Put the repository root and pipeline/ on the path before those
# imports run, so this script works when invoked directly as documented rather than
# only under pytest, whose conftest.py does the same thing.
_ROOT = _Path(__file__).resolve().parent.parent
for _entry in (_ROOT, _ROOT / "pipeline"):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

import argparse
import logging
import pandas as pd
import re
import os
import helpers as hp
from joblib import dump, load
from pathlib import Path
from dotenv import load_dotenv
import subprocess
import runpy
import shutil
import gzip
from typing import Tuple, Optional
import requests
import time
from Bio.Seq import Seq
from pyfaidx import Fasta
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from functools import partial
import json
import seaborn as sns
import matplotlib.pyplot as plt
import math
from itertools import combinations
import tempfile
from multiprocessing import Pool
from datetime import datetime
import numpy as np


# ── Phase selection ──
# Phases are sequential: preprocess → db → guides
# Use --start-from to skip earlier phases that already completed.
# Command example: python scripts/run_pipeline.py --start_from db
PHASE_ORDER = ["preprocess", "db", "guides", "plots"]
parser = argparse.ArgumentParser(description="ADAR/ASD pipeline")
parser.add_argument(
    "--start_from",
    choices=PHASE_ORDER,
    default="preprocess",
    help="Phase to start from (default: preprocess). Runs this phase and all subsequent ones.",
)
parser.add_argument(
    "--env-file",
    dest="env_file",
    default=os.environ.get("ASD_ENV_FILE", ".env"),
    help="Environment file naming the database to write to (default: .env, or $ASD_ENV_FILE). "
         "This selects the target database, so state it deliberately when a run must not "
         "touch the primary one.",
)
args = parser.parse_args()
START_INDEX = PHASE_ORDER.index(args.start_from)

PROJECT_ROOT = Path(__file__).resolve().parent.parent

def resolve(p):
    p = Path(p).expanduser()
    return p if p.is_absolute() else (PROJECT_ROOT / p).resolve()

# ── Load the environment ──
# The env file names the database this run writes to, so it is an explicit argument
# rather than a fixed path: a verification run must be able to target a scratch stack
# without editing this file.
env_path = resolve(args.env_file)
if not env_path.exists():
    raise RuntimeError(f"env file not found: {env_path}")
load_dotenv(env_path, override=True)

GENOME_REFERENCE_HG19     = resolve(os.environ["GENOME_REFERENCE_HG19"])
VARICARTA_PATH = resolve(os.environ["VARICARTA_PATH"])
OUTPUT_DIR     = resolve(os.environ["OUTPUT_DIR"])
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR      = OUTPUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
BLAT_RESULTS_DIR = OUTPUT_DIR / "blat_results"
BLAT_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
VEP_CACHE      = resolve(os.environ["VEP_CACHE"])
DB_DIR         = resolve(os.environ["DB_DIR"])
POSTGRES_DIR = resolve(os.environ["POSTGRES_DIR"])
POSTGRES_DIR.mkdir(parents=True, exist_ok=True)
POSTGRES_DB = os.environ["POSTGRES_DB"]
POSTGRES_USER = os.environ["POSTGRES_USER"]
POSTGRES_PASSWORD = os.environ["POSTGRES_PASSWORD"]
POSTGRES_PORT = os.environ["POSTGRES_PORT"]
POSTGRES_DOCKER_CONTAINER = os.environ["POSTGRES_DOCKER_CONTAINER"]
HASURA_URL = os.environ["HASURA_URL"]
HASURA_ADMIN_SECRET = os.environ["HASURA_ADMIN_SECRET"]
HASURA_DOCKER_CONTAINER = os.environ["HASURA_DOCKER_CONTAINER"]


# ── Logging Setup ──
LOG_DIR = OUTPUT_DIR / "logs"
os.makedirs(LOG_DIR, exist_ok=True)
PREPROCESS_LOG_FILE = LOG_DIR / "preprocess.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(PREPROCESS_LOG_FILE, mode="w", encoding="utf-8"),
    ],
)

def switch_log_file(log_path):
    """Replace the FileHandler on the root logger to write to a new file."""
    logger = logging.getLogger()
    for handler in logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            handler.close()
            logger.removeHandler(handler)
    new_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    new_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    logger.addHandler(new_handler)
    logging.info(f"Logging to: {log_path}")

logging.info(f"Loaded .env from: {env_path}")

HASURA_URL = str(os.environ.get("HASURA_URL"))
HASURA_ADMIN_SECRET = str(os.environ.get("HASURA_ADMIN_SECRET"))

HEADERS = {
    "x-hasura-admin-secret": HASURA_ADMIN_SECRET
}

gql = hp.gql

logging.info(f"Pipeline starting from phase: {args.start_from}")

# ── Shared paths (needed by multiple phases) ──
path_varicarta_standard_positive = OUTPUT_DIR / "varicarta_standard_positive_clean_variants.vcf"
vep_output_file_name = "varicarta_vepped"
vep_output_dir = OUTPUT_DIR / "VEP"
vep_output_dir.mkdir(parents=True, exist_ok=True)

# --- the contract between the runner and a stage -----------------------------
# Stages used to be exec()'d, which handed each one this module's entire scope:
# whatever a stage happened to use, it silently got. What it actually needed was
# invisible, so no stage could be read or run on its own.
#
# runpy.run_path with explicit init_globals inverts that. The names below are the
# whole interface -- derived by static analysis of what each stage reads and never
# defines -- and anything a stage needs beyond them now fails with a NameError
# naming the missing value, instead of quietly resolving to whatever the runner
# happened to have in scope.
stage_globals = {
    "hp": hp,
    "gql": gql,
    "PROJECT_ROOT": PROJECT_ROOT,
    "OUTPUT_DIR": OUTPUT_DIR,
    "DB_DIR": DB_DIR,
    "BLAT_RESULTS_DIR": BLAT_RESULTS_DIR,
    "VARICARTA_PATH": VARICARTA_PATH,
    "GENOME_REFERENCE_HG19": GENOME_REFERENCE_HG19,
    "VEP_CACHE": VEP_CACHE,
    "HASURA_URL": HASURA_URL,
    "HASURA_ADMIN_SECRET": HASURA_ADMIN_SECRET,
    "POSTGRES_DB": POSTGRES_DB,
    "POSTGRES_USER": POSTGRES_USER,
    "POSTGRES_DOCKER_CONTAINER": POSTGRES_DOCKER_CONTAINER,
    "resolve": resolve,
    "path_varicarta_standard_positive": path_varicarta_standard_positive,
    "vep_output_dir": vep_output_dir,
    "vep_output_file_name": vep_output_file_name,
}


def run_stage(filename, label, log_file=None):
    """Run one pipeline stage with exactly the shared state it declares."""
    if log_file is not None:
        switch_log_file(log_file)
    logging.info(f"=== Running {label} ===")
    runpy.run_path(str(PROJECT_ROOT / "pipeline" / filename), init_globals=dict(stage_globals))


if START_INDEX <= 0:
    run_stage("phase1_preprocess.py", "Phase 1: Preprocess")
else:
    logging.info("Skipping Phase 1: Preprocess")

if START_INDEX <= 1:
    run_stage("phase2_db.py", "Phase 2: DB Insertion", LOG_DIR / "db_insertion.log")
else:
    logging.info("Skipping Phase 2: DB Insertion")

if START_INDEX <= 2:
    run_stage("phase3_guides.py", "Phase 3: Guides", LOG_DIR / "guides.log")
else:
    logging.info("Skipping Phase 3: Guides")


logging.info("Pipeline finished.")
