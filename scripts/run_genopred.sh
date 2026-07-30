#!/usr/bin/env bash
# Purpose: prepare and execute the eight-trait GenoPred/PRS-CS-auto analysis.
# Inputs: project configuration, standardised GWAS files, validated target pfiles,
#         custom GenoPred reference and fixed GenoPred runtime.
# Outputs: GenoPred pipeline directory, logs and a compact validation record.
# Usage: bash scripts/run_genopred.sh --config config/project.yml [--dry-run]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG=""
CORES="${SLURM_CPUS_PER_TASK:-4}"
DRY_RUN=0
PREPARE_ONLY=0
SKIP_COMMIT_CHECK=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --skip-commit-check) SKIP_COMMIT_CHECK=1; shift ;;
    -h|--help)
      sed -n '1,12p' "$0"
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "--config is required" >&2
  exit 2
fi
if ! [[ "$CORES" =~ ^[1-9][0-9]*$ ]]; then
  echo "--cores must be a positive integer" >&2
  exit 2
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON_BIN" -m gendep.prs.genopred \
  --project-config "$CONFIG" \
  --genopred-config "$ROOT/config/genopred.yml" \
  --gwas-config "$ROOT/config/gwas.yml" \
  --traits-file "$ROOT/config/traits.tsv"

RUN_ENV="$($PYTHON_BIN - "$CONFIG" <<'PY'
import sys
from pathlib import Path
import yaml
config=Path(sys.argv[1]).resolve()
d=yaml.safe_load(config.read_text())
root=Path(d.get('project_root', config.parent.parent)).expanduser()
if not root.is_absolute():
    root=(config.parent/root).resolve()
out=Path(d['outputs']['genopred_run']).expanduser()
if not out.is_absolute():
    out=root/out
print((out/'config'/'run.env').resolve())
PY
)"
# The generated file contains only shell-quoted paths and fixed commit strings.
# shellcheck disable=SC1090
source "$RUN_ENV"

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/manifests"
LOG="$RUN_ROOT/logs/genopred_workflow.log"
VALIDATION="$RUN_ROOT/manifests/genopred_execution_validation.txt"

for executable in "$SNAKEMAKE_BIN" "$RSCRIPT_BIN" "$RUNTIME_PYTHON" "$PLINK_BIN" "$PLINK2_BIN"; do
  [[ -x "$executable" ]] || { echo "Executable is missing: $executable" >&2; exit 1; }
done
[[ -s "$RUN_CONFIG" ]] || { echo "Run configuration is missing: $RUN_CONFIG" >&2; exit 1; }

if [[ "$SKIP_COMMIT_CHECK" -eq 0 ]]; then
  observed_commit="$(git -C "$GENOPRED_ROOT" rev-parse HEAD)"
  [[ "$observed_commit" == "$EXPECTED_GENOPRED_COMMIT" ]] || {
    echo "GenoPred commit mismatch: $observed_commit" >&2
    exit 1
  }
  observed_genoutils="$($RSCRIPT_BIN -e 'd <- packageDescription("GenoUtils"); cat(d$RemoteSha)')"
  [[ "$observed_genoutils" == "$EXPECTED_GENOUTILS_COMMIT" ]] || {
    echo "GenoUtils commit mismatch: $observed_genoutils" >&2
    exit 1
  }
else
  observed_commit="NOT_CHECKED"
  observed_genoutils="NOT_CHECKED"
fi

if [[ "$PREPARE_ONLY" -eq 1 ]]; then
  echo "GENOPRED_PREPARATION=PASS"
  exit 0
fi

RUNTIME_BIN="$(dirname "$RUNTIME_PYTHON")"
SNAKEMAKE_BIN_DIR="$(dirname "$SNAKEMAKE_BIN")"
export PATH="$RUNTIME_BIN:$SNAKEMAKE_BIN_DIR:/usr/bin:/bin"
export CONDA_DEFAULT_ENV="genopred"
export CONDA_PREFIX="$(dirname "$RUNTIME_BIN")"
export PYTHONNOUSERSITE=1
export R_LIBS_USER=""
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
hash -r

cd "$GENOPRED_PIPELINE"
command=("$SNAKEMAKE_BIN" --snakefile Snakefile --configfile "$RUN_CONFIG" --cores "$CORES" --printshellcmds --rerun-incomplete --latency-wait 120 --show-failed-logs)
if [[ "$DRY_RUN" -eq 1 ]]; then
  command+=(--dry-run)
fi
command+=(output_all)

set +e
"${command[@]}" >"$LOG" 2>&1
status=$?
set -e

postflight="NOT_RUN"
if [[ "$status" -eq 0 && "$DRY_RUN" -eq 0 ]]; then
  postflight="PASS"
  for trait in MDD ANX BIP SCZ NEUR INSOM SWB EA; do
    cleaned="$PIPELINE_OUTPUT/reference/gwas_sumstat/$trait/$trait-cleaned.gz"
    score="$PIPELINE_OUTPUT/reference/pgs_score_files/prscs/$trait/ref-$trait.score.gz"
    profile="$PIPELINE_OUTPUT/GENDEP/pgs/TRANS/prscs/$trait/GENDEP-$trait-TRANS.profiles"
    for required in "$cleaned" "$score" "$profile"; do
      if [[ ! -s "$required" ]]; then
        echo "Required GenoPred output is missing or empty: $required" >&2
        postflight="FAIL"
      fi
    done
    [[ -s "$cleaned" ]] && gzip -t "$cleaned"
    [[ -s "$score" ]] && gzip -t "$score"
  done
  [[ "$postflight" == "PASS" ]] || status=1
fi

cat > "$VALIDATION" <<EOF
run_config=$RUN_CONFIG
pipeline_output=$PIPELINE_OUTPUT
cores=$CORES
dry_run=$DRY_RUN
genopred_commit=$observed_commit
genoutils_commit=$observed_genoutils
workflow_exit_status=$status
postflight=$postflight
workflow_log=$LOG
validation=$([[ "$status" -eq 0 ]] && echo PASS || echo FAIL)
EOF

cat "$VALIDATION"
[[ "$status" -eq 0 ]] || exit "$status"
echo "GENOPRED_EXECUTION=PASS"
