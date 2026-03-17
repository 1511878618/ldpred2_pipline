#!/usr/bin/env bash
set -euo pipefail

########################################
# Defaults (optional hyperparameters)
########################################
LDCache_default="LDCache/chr#."
pfile_default="/pmaster/xutingfeng/dataset/ukb/dataset/snp/ukb_imputed_v3_qc/qc_pgen_hg19_hapmap3/ukb_imp_v3_hg37_qc_hapmap3_db150"
tuning_set_pheno_col_default="3"
score_col_nums_default="4-101"
score_col_start_default="5"

########################################
# Parse CLI
########################################
usage() {
  cat <<'EOF'
Usage:
  bash run_prs.sh -g GWAS -f TUNING -o OUTDIR [options]

Required:
  -g <file>   GWAS file (ldpred2 input)
  -f <file>   tuning set regenie file
  -o <dir>    output directory

Optional:
  -L <path>   LDCache pattern (default: LDCache/chr#.)
  -p <pfile>  plink2 --pfile prefix
  -t <int>    tuning set phenotype column (default: 3)

  --score-col-nums   <range>  (default: 4-101)
  --score-col-start  <int>    (default: 5)

  --dry-run
  -h, --help

Example:
  bash run_prs.sh \
    -g Input/MVP/GWAS/xxx.tsv.gz \
    -f Input/MVP/TuningSet/mi.regenie \
    -t 3 \
    -o Output/MVP/mi
EOF
}

# initialize as EMPTY → 强制用户输入
gwas_file=""
tuning_set_file=""
output_dir=""

# optional with defaults
LDCache="$LDCache_default"
pfile="$pfile_default"
tuning_set_pheno_col="$tuning_set_pheno_col_default"
score_col_nums="$score_col_nums_default"
score_col_start="$score_col_start_default"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) gwas_file="$2"; shift 2 ;;
    -f) tuning_set_file="$2"; shift 2 ;;
    -o) output_dir="$2"; shift 2 ;;
    -L) LDCache="$2"; shift 2 ;;
    -p) pfile="$2"; shift 2 ;;
    -t) tuning_set_pheno_col="$2"; shift 2 ;;

    --score-col-nums) score_col_nums="$2"; shift 2 ;;
    --score-col-start) score_col_start="$2"; shift 2 ;;

    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

########################################
# Mandatory argument checks
########################################
[[ -n "$gwas_file" ]] || { echo "ERROR: -g <gwas_file> is required" >&2; usage; exit 1; }
[[ -n "$tuning_set_file" ]] || { echo "ERROR: -f <tuning_set_file> is required" >&2; usage; exit 1; }
[[ -n "$output_dir" ]] || { echo "ERROR: -o <output_dir> is required" >&2; usage; exit 1; }

mkdir -p "$output_dir"

########################################
# Helpers
########################################
log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

need_file() { [[ -s "$1" ]] || die "Missing or empty input: $1"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }

run_step() {
  local name="$1"; shift
  local out="$1"; shift

  if [[ -s "$out" ]]; then
    log "[SKIP] $name (exists: $out)"
    return 0
  fi

  log "[RUN ] $name"
  if [[ "$dry_run" -eq 1 ]]; then
    printf '  '; printf '%q ' "$@"; echo
    return 0
  fi

  "$@"
  [[ -s "$out" ]] || die "$name finished but output missing/empty: $out"
  log "[DONE] $name -> $out"
}




########################################
# Preflight
########################################
need_file "$gwas_file"
need_file "$tuning_set_file"

need_cmd Rscript
need_cmd plink2
need_cmd python
need_cmd zcat
need_cmd gzip
need_cmd csvcut
need_cmd csvformat

startTime=$(date +%s)
output_dir="${output_dir%%/}/"

log "Run config:"
log "  gwas_file            = $gwas_file"
log "  tuning_set_file      = $tuning_set_file"
log "  tuning_set_pheno_col = $tuning_set_pheno_col"
log "  output_dir           = $output_dir"

########################################
# Steps (保持你原有逻辑)
########################################
step1_out="${output_dir}/ldpred2_grid.tsv.gz"

run_step "LDpred2 weights" "$step1_out" \
  Rscript LDpred2AndLassosum2LDCaches_v1.R -o "$output_dir" -g "$gwas_file" -p "$LDCache"

# Run step2; assign rsid
step2_out="${output_dir}/ldpred2_grid_rsid.tsv.gz"
run_step "Assign rsid" "$step2_out" \
  bash -c "
    zcat '$step1_out' \
      | python KeyMapReplacer_v2.py -p EUR_38.rsid -k 1 -r 1 \
      | gzip -c > '$step2_out'
  "

# -----------------------------
# Step 2.5: NA check + replace to 0 (only weight columns)
# -----------------------------
step2_clean_out="${output_dir}/ldpred2_grid_rsid.na0.tsv.gz"

run_step "NA check + replace weights to 0" "$step2_clean_out" \
  bash -c "
    zcat \"$step2_out\" | \
    awk -F'\t' -v OFS='\t' '{
      for (i=4; i<=NF; i++)
        if (\$i==\"\" || \$i==\"NA\" || \$i==\"NaN\" || \$i==\".\") \$i=0
      print
    }' | gzip > \"$step2_clean_out\"
  "

# step3
step3_out="${output_dir}/LDPred2_PRS.sscore.gz"

run_step "plink2 score + gzip" "$step3_out" \
  bash -c "
    plink2 \
      --score '$step2_clean_out' 1 2 header-read cols=+scoresums \
      --score-col-nums '$score_col_nums' \
      --pfile '$pfile' \
      --out '${output_dir}/LDPred2_PRS' \
    && gzip -f '${output_dir}/LDPred2_PRS.sscore'
  "


step4_out="${output_dir}/LDPred2_PRS.auc.tsv"
run_step "Evaluate AUC" "$step4_out" \
  python cal_auc.py \
    -f "$tuning_set_file" \
    -p "$tuning_set_pheno_col" \
    -s "$step3_out" \
    --score-col-start "$score_col_start" \
    -o "$step4_out"

########################################
# Step 5: pick best score + extract best column (from step1_out + step4_out)
########################################
step5_out="${output_dir}/ldpred2_grid.best.tsv.gz"

if [[ -s "$step5_out" ]]; then
  log "[SKIP] Extract best weights (exists: $step5_out)"
else
  log "[RUN ] Extract best weights"

  # inputs check
  need_file "$step4_out"
  need_file "$step2_clean_out"

  # parse best score name from auc file (2nd row, 2nd column)
  best_score_name="$(
    awk -F'\t' 'NR==2 {print $2}' "$step4_out" \
      | sed 's/_SUM$//; s/_AVG$//'
  )"
  [[ -n "$best_score_name" ]] || die "Failed to parse best_score_name from: $step4_out"
  log "Best score column: $best_score_name"

  if [[ "$dry_run" -eq 1 ]]; then
    log "  (dry-run) would run: zcat '$step1_out' | csvcut -t -c 1,2,3,'$best_score_name' | csvformat -T | gzip -c > '$step5_out'"
  else
    zcat "$step2_clean_out" \
      | csvcut -t -c 1,2,"$best_score_name" \
      | csvformat -T \
      | gzip -c > "$step5_out"

    [[ -s "$step5_out" ]] || die "Extract best weights finished but output missing/empty: $step5_out"
    log "[DONE] Extract best weights -> $step5_out"
  fi
fi



########################################
# Time log
########################################
EndTime=$(date +%s)
if [[ "$dry_run" -eq 0 ]]; then
  {
    echo "Start: $(date -d "@$startTime" '+%F %T')"
    echo "End:   $(date -d "@$EndTime" '+%F %T')"
    echo "TotalSeconds: $((EndTime - startTime))"
  } > "${output_dir}/time.log"
fi

log "All done."
