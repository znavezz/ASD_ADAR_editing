#!/usr/bin/env bash
# reset_db.sh — Revert the workspace to the state it was in before a given pipeline phase.
#
# Usage:
#   ./reset_db.sh [--before <1|2|3>] [--env-file <path>]
#
#   --before 1  Full reset: remove all pipeline output (default when no args given)
#   --before 2  Remove DB data, staging CSVs, and all guide output (keep phase 1 outputs)
#   --before 3  Remove only guide/bystander rows from the DB and related files
#
#   --env-file  Environment file naming the database to reset. Defaults to
#               .env.test (or $ASD_ENV_FILE): this script DESTROYS data, so the
#               scratch stack is the default and the real one must be asked for.
#               Resetting the published database additionally requires
#               ASD_ALLOW_PUBLISHED_RESET=yes.
#
# Each option removes everything produced by the specified phase and all later phases.
set -euo pipefail

cd "$(dirname "$0")"

# Default to the scratch stack. This script deletes data, and the alternative default
# was the database that produced the manuscript: `./reset_db.sh` with no arguments
# destroyed it. Choosing the destructive target has to be a deliberate act, so the
# safe stack is what you get for free.
ENV_FILE="${ASD_ENV_FILE:-.env.test}"
_args=(); while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2 ;;
    *) _args+=("$1"); shift ;;
  esac
done
set -- ${_args+"${_args[@]}"}

[[ -f "$ENV_FILE" ]] || { echo "env file not found: $ENV_FILE" >&2; exit 1; }

# Export, don't just source. `docker compose` does not read the shell's unexported
# variables: it reads ./.env itself. Sourcing .env.test without exporting therefore
# truncated the test database (docker exec uses the shell variable) while resolving
# COMPOSE_PROJECT_NAME to the production project, so `docker compose down` would have
# stopped the wrong stack.
set -a
source "$ENV_FILE"
set +a

echo "Target: ${POSTGRES_DOCKER_CONTAINER} / ${POSTGRES_DB} (${ENV_FILE})"
echo "Compose project: $(docker compose config 2>/dev/null | awk '/^name:/{print $2; exit}')"

# The published database is the one that produced the manuscript. A full reset deletes
# its data directory outright, and --before 1 is the default when no arguments are given,
# so refuse unless the caller says so in as many words. ASD_PUBLISHED_CONTAINER names the
# instance to protect; set it empty to disable the guard entirely.
: "${ASD_PUBLISHED_CONTAINER:=asd_adar_postgres}"
if [[ -n "${ASD_PUBLISHED_CONTAINER}" \
   && "${POSTGRES_DOCKER_CONTAINER}" == "${ASD_PUBLISHED_CONTAINER}" \
   && "${ASD_ALLOW_PUBLISHED_RESET:-}" != "yes" ]]; then
  echo "" >&2
  echo "Refusing: ${POSTGRES_DOCKER_CONTAINER} is the published database." >&2
  echo "  A reset deletes ${POSTGRES_DIR}, which is the data behind the paper." >&2
  echo "" >&2
  echo "  To reset a scratch stack instead:  $0 --before 1 --env-file .env.test" >&2
  echo "  To proceed anyway:                 ASD_ALLOW_PUBLISHED_RESET=yes $0 $* --env-file ${ENV_FILE}" >&2
  exit 1
fi

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
  echo "Usage: $0 [--before <1|2|3>] [--env-file <path>]"
  echo ""
  echo "  --before 1  Full reset — remove all pipeline output (default)"
  echo "  --before 2  Remove DB, staging CSVs, and guide output (keep VCF/VEP)"
  echo "  --before 3  Remove only guide/bystander data from DB and related files"
  exit 1
}

BEFORE=1  # default: full reset (backward-compatible with original behaviour)

if [[ $# -gt 0 ]]; then
  if [[ "$1" == "--before" ]]; then
    [[ $# -lt 2 ]] && { echo "Error: --before requires a phase number (1, 2, or 3)" >&2; usage; }
    BEFORE="$2"
    shift 2
  else
    echo "Error: unknown argument: $1" >&2
    usage
  fi
fi

[[ "$BEFORE" == "1" || "$BEFORE" == "2" || "$BEFORE" == "3" ]] \
  || { echo "Error: phase must be 1, 2, or 3 (got: ${BEFORE})" >&2; usage; }

[[ $# -eq 0 ]] || { echo "Error: unexpected arguments: $*" >&2; usage; }

# ── Helper: ensure containers are running (needed to TRUNCATE) ────────────────
ensure_db_running() {
  local running
  running=$(docker ps --filter "name=^${POSTGRES_DOCKER_CONTAINER}$" --format '{{.Names}}')
  if [[ -z "$running" ]]; then
    echo "Starting Docker services..."
    USER_ID=$(id -u) GROUP_ID=$(id -g) docker compose up -d
    echo "Waiting for Postgres to be ready..."
    local i=0
    while ! docker exec "${POSTGRES_DOCKER_CONTAINER}" \
        pg_isready -U "${POSTGRES_USER}" -d "${POSTGRES_DB}" -q 2>/dev/null; do
      # NOT ((i++)): with i=0 that evaluates to 0, which bash reports as exit status 1,
      # and `set -e` then kills the script on the first iteration. This loop never ran
      # more than once, so the branch that starts a stopped container never worked.
      i=$((i + 1))
      [[ $i -ge 30 ]] && { echo "Error: Postgres did not become ready in time." >&2; exit 1; }
      sleep 2
    done
    echo "  Postgres is ready."
  fi
}

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 3 cleanup — guide/bystander rows and BLAT/plot files
# (always runs, since phase 3 is the last phase)
# ═══════════════════════════════════════════════════════════════════════════════
echo "=== Reverting Phase 3 outputs ==="

ensure_db_running

echo "Truncating guide and bystander tables..."
docker exec "${POSTGRES_DOCKER_CONTAINER}" \
  psql -U "${POSTGRES_USER}" -d "${POSTGRES_DB}" \
  -c "TRUNCATE guides CASCADE;" \
  -c "TRUNCATE neighbor_edits CASCADE;"
echo "  ✓ guides, variants_guides, bystanders, bystanders_consequences, neighbor_edits (+ scores/consequences) cleared"

echo "Removing BLAT results..."
rm -rf "${OUTPUT_DIR}/blat_results/"*

echo "Removing plots..."
rm -rf "${OUTPUT_DIR}/plots/"*

rm -f "${OUTPUT_DIR}/logs/guides.log"
echo "  ✓ Phase 3 file outputs removed"

if [[ "${BEFORE}" == "3" ]]; then
  echo ""
  echo "Done. Re-run from phase 3:  python scripts/run_pipeline.py --start_from guides"
  exit 0
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 2 cleanup — Docker containers, Postgres data, staging CSVs
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Reverting Phase 2 outputs ==="

echo "Stopping Docker services..."
USER_ID=$(id -u) GROUP_ID=$(id -g) docker compose down

echo "Removing Postgres data directory..."
rm -rf "${POSTGRES_DIR:?}"/*

echo "Removing staging CSV files..."
rm -f "${DB_DIR}/staging_varicarta_variants_raw.csv"
rm -f "${DB_DIR}/staging_vep_annotations_raw.csv"

rm -f "${OUTPUT_DIR}/logs/db_insertion.log"
echo "  ✓ Phase 2 outputs removed"

if [[ "${BEFORE}" == "2" ]]; then
  echo ""
  echo "Done. Re-run from phase 2:  python scripts/run_pipeline.py --start_from db"
  exit 0
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1 cleanup — preprocessed VCF and VEP output
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Reverting Phase 1 outputs ==="

echo "Removing preprocessed VCF..."
rm -f "${OUTPUT_DIR}/varicarta_standard_positive_clean_variants.vcf"

echo "Removing VEP output directory..."
rm -rf "${OUTPUT_DIR}/VEP"

rm -f "${OUTPUT_DIR}/logs/preprocess.log"
echo "  ✓ Phase 1 outputs removed"

echo ""
echo "Done. Re-run from the beginning:  python scripts/run_pipeline.py"
