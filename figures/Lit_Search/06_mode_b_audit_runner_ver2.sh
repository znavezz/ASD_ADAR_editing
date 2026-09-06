#!/usr/bin/env bash
# =============================================================================
# Ver-2 Mode B audit runner — Sonnet extract "which of the 631 variants are
# named in this paper" per (gene, PMID) prompt.  Concurrency 5, resumable.
#
# Launch inside `screen -dmS mode_b_v2 bash /path/to/this.sh` so it survives.
# =============================================================================
set -uo pipefail

VER2="${ASD_SOURCE}"
PROMPTS_DIR="$VER2/Intermediate_tables/Lit_Search/mode_b_ver2/prompts"
OUT_DIR="$VER2/Results_Plots_Tables/Lit_Search/mode_b_ver2/audits"
LOG_DIR="$VER2/Intermediate_tables/Lit_Search/mode_b_ver2/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

CONCURRENCY=5
MODEL="claude-sonnet-4-6"
MAX_BUDGET=0.5

RUN_LOG="$LOG_DIR/run_$(date +%Y%m%d_%H%M%S).log"
echo "=== ver-2 mode B audit start $(date) ===" | tee -a "$RUN_LOG"
NAMES=$(ls "$PROMPTS_DIR" 2>/dev/null | sed 's/\.prompt$//')
TOTAL=$(echo "$NAMES" | wc -l)
echo "Total prompts: $TOTAL  Concurrency: $CONCURRENCY  Budget/call: \$$MAX_BUDGET" | tee -a "$RUN_LOG"

process_one() {
  local name="$1"
  local prompt="$PROMPTS_DIR/$name.prompt"
  local out="$OUT_DIR/$name.txt"
  local err="$LOG_DIR/$name.err"
  if [[ -f "$out" && -s "$out" ]]; then echo "[skip] $name"; return 0; fi
  local t0=$(date +%s)
  claude --print --no-session-persistence --permission-mode bypassPermissions \
         --model "$MODEL" --max-budget-usd "$MAX_BUDGET" --output-format text \
         < "$prompt" > "$out" 2> "$err"
  local rc=$?; local t1=$(date +%s); local dur=$((t1 - t0))
  if [[ $rc -eq 0 && -s "$out" ]]; then
    echo "[done] $name (${dur}s)"; rm -f "$err"
  else
    echo "[FAIL] $name (rc=$rc, ${dur}s) — see $err"
    [[ -f "$out" ]] && mv "$out" "$out.failed.$(date +%s)"
    return 1
  fi
}
export -f process_one
export PROMPTS_DIR OUT_DIR LOG_DIR MODEL MAX_BUDGET

echo "$NAMES" | xargs -P "$CONCURRENCY" -I {} bash -c 'process_one "$@"' _ {} | tee -a "$RUN_LOG"

DONE=$(find "$OUT_DIR" -maxdepth 1 -type f -name '*.txt' ! -name '*.failed.*' 2>/dev/null | wc -l)
FAIL=$(find "$OUT_DIR" -maxdepth 1 -type f -name '*.failed.*' 2>/dev/null | wc -l)
echo "=== finished $(date) — done $DONE/$TOTAL  failed $FAIL ===" | tee -a "$RUN_LOG"
