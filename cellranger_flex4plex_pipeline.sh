#!/usr/bin/env bash
# =============================================================================
# CellRanger pipeline — GEM-X Flex Gene Expression, 4-plex (probe-based)
#
# Kit used: GEM-X Flex Gene Expression Human 4-plex (PN-1000793)
# Protocol:  CG000787 Rev B
# CellRanger: requires >= 7.0 (>= 10.0 recommended for GEM-X Flex v2)
#
# How GEM-X Flex multiplexing works (differs from CMO/Hashtag):
#   - Each sample is hybridised with a unique Probe Barcode (BC001–BC004)
#   - All 4 samples are pooled into ONE GEM well → ONE pair of FASTQ files
#   - CellRanger demultiplexes by Probe Barcode, NOT by a separate CMO library
#   - There is NO separate multiplexing-capture FASTQ; only Gene Expression FASTQs
#   - The correct pipeline is: cellranger multi with a [samples] section
#
# FASTQ files provided: SO_14534_R1.fastq.gz / SO_14534_R2.fastq.gz
#
# Usage:
#   chmod +x cellranger_flex4plex_pipeline.sh
#   ./cellranger_flex4plex_pipeline.sh [--dry-run] [options]
#
# =============================================================================

set -euo pipefail

# ─── USER CONFIGURATION ──────────────────────────────────────────────────────
# Edit these before running.

# CellRanger genome reference
TRANSCRIPTOME="/home/cxbl/hari/snrna_project/reference/refdata-gex-GRCh38-2024-A"

# Probe set CSV bundled with your CellRanger installation.
# IMPORTANT: Based on the NGS Library Report (CG000787 Rev B protocol),
# this was prepared with Chromium Fixed RNA Profiling (Flex v1) chemistry.
# Confirm with Genotypic which probe set version was used.
# For Flex v1 use:
#   cellranger-X.Y.Z/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv
# For Flex v1 use:
#   cellranger-X.Y.Z/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv
PROBE_SET="/home/cxbl/hari/snrna_project/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"

# Directory containing SO_14534_R1.fastq.gz and SO_14534_R2.fastq.gz
FASTQ_DIR="/home/cxbl/hari/snrna_project/raw_reads"

# The FASTQ prefix (the part before _R1/_R2) — used as fastq_id in config
FASTQ_ID="SO_14534"  # Matches SO_14534_R1.fastq.gz / SO_14534_R2.fastq.gz

# Output directory
OUTPUT_DIR="/home/cxbl/hari/snrna_project/cellranger_out"

# Run ID (name for the cellranger multi output folder)
RUN_ID="SO_14534_flex4plex"

# Resources
CORES=16
MEM_GB=128

# ─── SAMPLE CONFIGURATION ────────────────────────────────────────────────────
# Map each sample name to its Probe Barcode ID (BC001–BC004 for the 4-plex kit)
# BC001 = Human WTA Probes BC001 (PN-2001259)
# BC002 = Human WTA Probes BC002 (PN-2001260)
# BC003 = Human WTA Probes BC003 (PN-2001261)
# BC004 = Human WTA Probes BC004 (PN-2001262)
#
# The mapping of sample→barcode must match Table 2 of your NGS Library Report.
# Adjust sample names (SC001–SC004) to your actual sample identifiers.

declare -A SAMPLE_TO_BC=(
  ["SC001"]="BC001"   # SO_14534_SC001 → Human WTA Probe BC001
  ["SC004"]="BC002"   # SO_14534_SC004 → Human WTA Probe BC002
  ["SC006"]="BC003"   # SO_14534_SC006 → Human WTA Probe BC003
  ["SC008"]="BC004"   # SO_14534_SC008 → Human WTA Probe BC004
)

# Ordered list (bash associative arrays don't preserve order)
SAMPLE_ORDER=("SC001" "SC004" "SC006" "SC008")

# ─── FLAGS ───────────────────────────────────────────────────────────────────
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)       DRY_RUN=true;          shift ;;
    --transcriptome) TRANSCRIPTOME="$2";    shift 2 ;;
    --probe-set)     PROBE_SET="$2";        shift 2 ;;
    --fastq-dir)     FASTQ_DIR="$2";        shift 2 ;;
    --fastq-id)      FASTQ_ID="$2";         shift 2 ;;
    --output-dir)    OUTPUT_DIR="$2";       shift 2 ;;
    --run-id)        RUN_ID="$2";           shift 2 ;;
    --cores)         CORES="$2";            shift 2 ;;
    --mem)           MEM_GB="$2";           shift 2 ;;
    -h|--help)
      grep '^#' "$0" | head -30 | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ─── HELPERS ─────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die()  { echo "[ERROR] $*" >&2; exit 1; }
run()  {
  if $DRY_RUN; then
    echo "[DRY-RUN] $*"
  else
    "$@"
  fi
}

# ─── PRE-FLIGHT CHECKS ───────────────────────────────────────────────────────
preflight() {
  log "=== Pre-flight checks ==="

  command -v cellranger >/dev/null 2>&1 || die "cellranger not found on PATH"
  CR_VER=$(cellranger --version 2>&1 | grep -oP '\d+\.\d+\.\d+' | head -1 || echo "unknown")
  log "CellRanger version: ${CR_VER}"

  # GEM-X Flex v2 requires CR >= 10.0; Flex v1 works from CR 7.0
  CR_MAJOR=$(echo "$CR_VER" | cut -d. -f1)
  if [[ "$CR_MAJOR" -lt 7 ]]; then
    die "CellRanger >= 7.0 required for Flex; found ${CR_VER}"
  fi
  if [[ "$CR_MAJOR" -lt 10 ]]; then
    log "[WARN] CellRanger >= 10.0 recommended for GEM-X Flex (v2). Found ${CR_VER}."
    log "       If your data is Flex v2, upgrade to CR 10+."
  fi

  [[ -d "$TRANSCRIPTOME" ]] || die "Transcriptome not found: ${TRANSCRIPTOME}"
  [[ -f "$PROBE_SET" ]]     || die "Probe set CSV not found: ${PROBE_SET}"
  [[ -d "$FASTQ_DIR" ]]     || die "FASTQ directory not found: ${FASTQ_DIR}"

  # Check R1/R2 files are present
  local R1="${FASTQ_DIR}/${FASTQ_ID}_R1.fastq.gz"
  local R2="${FASTQ_DIR}/${FASTQ_ID}_R2.fastq.gz"
  [[ -f "$R1" ]] || die "R1 FASTQ not found: ${R1}"
  [[ -f "$R2" ]] || die "R2 FASTQ not found: ${R2}"
  log "FASTQ files found: ${R1##*/}  ${R2##*/}"

  mkdir -p "$OUTPUT_DIR"
  log "Output directory: ${OUTPUT_DIR}"
  log "Pre-flight checks passed."
}

# ─── GENERATE multi config.csv ───────────────────────────────────────────────
write_config() {
  local CONFIG="${OUTPUT_DIR}/flex4plex_config.csv"
  log "Writing cellranger multi config: ${CONFIG}"

  cat > "$CONFIG" <<CSV
[gene-expression]
reference,${TRANSCRIPTOME}
probe-set,${PROBE_SET}
create-bam,false
# expect-cells: set to total expected cells across all 4 samples combined
# 4 samples × 80,000 targeted nuclei = 320,000 total
# expect-cells,320000

[libraries]
fastq_id,fastqs,feature_types
${FASTQ_ID},${FASTQ_DIR},Gene Expression

[samples]
sample_id,probe_barcode_ids
CSV

  for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    echo "${SAMPLE},${SAMPLE_TO_BC[$SAMPLE]}" >> "$CONFIG"
  done

  log "Config written:"
  cat "$CONFIG"
  echo "$CONFIG"
}

# ─── RUN cellranger multi ────────────────────────────────────────────────────
run_cellranger_multi() {
  local CONFIG="$1"

  log "=== Running cellranger multi ==="
  log "Run ID:   ${RUN_ID}"
  log "Samples:  ${SAMPLE_ORDER[*]}"
  log "Barcodes: ${SAMPLE_TO_BC[*]}"

  # Check if already completed
  local DONE_FLAG="${OUTPUT_DIR}/${RUN_ID}/outs/per_sample_outs/${SAMPLE_ORDER[-1]}/metrics_summary.csv"
  if [[ -f "$DONE_FLAG" ]]; then
    log "[SKIP] cellranger multi output already exists. Delete ${OUTPUT_DIR}/${RUN_ID} to rerun."
    return 0
  fi

  cd "$OUTPUT_DIR"
  run cellranger multi \
    --id          "${RUN_ID}" \
    --csv         "${CONFIG}" \
    --localcores  "${CORES}" \
    --localmem    "${MEM_GB}"

  log "cellranger multi complete."
}

# ─── QC SUMMARY ──────────────────────────────────────────────────────────────
qc_summary() {
  log "=== QC summary ==="

  local PASS=0
  local FAIL=0

  printf "\n%-10s  %-12s  %-12s  %-10s  %-10s\n" \
    "Sample" "Est.Cells" "Median_Genes" "Pct_mito" "Saturation"
  printf "%-10s  %-12s  %-12s  %-10s  %-10s\n" \
    "------" "---------" "------------" "--------" "----------"

  for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    local METRICS="${OUTPUT_DIR}/${RUN_ID}/outs/per_sample_outs/${SAMPLE}/metrics_summary.csv"
    local WEB="${OUTPUT_DIR}/${RUN_ID}/outs/per_sample_outs/${SAMPLE}/web_summary.html"

    if [[ -f "$METRICS" ]]; then
      # Extract fields — column positions vary slightly across CR versions
      CELLS=$(awk -F',' 'NR==2{print $1}' "$METRICS" 2>/dev/null | tr -d '"' || echo "N/A")
      GENES=$(grep -i "median genes per cell" "$METRICS" | awk -F',' '{print $NF}' | tr -d '"\r' || echo "N/A")
      SAT=$(grep -i "sequencing saturation" "$METRICS" | awk -F',' '{print $NF}' | tr -d '"\r' || echo "N/A")
      MITO=$(grep -i "median.*mito" "$METRICS" | awk -F',' '{print $NF}' | tr -d '"\r' || echo "N/A")
      printf "%-10s  %-12s  %-12s  %-10s  %-10s\n" \
        "$SAMPLE" "$CELLS" "$GENES" "$MITO" "$SAT"
      PASS=$((PASS+1))
    else
      printf "%-10s  [metrics_summary.csv not found]\n" "$SAMPLE"
      FAIL=$((FAIL+1))
    fi

    [[ -f "$WEB" ]] && log "  Web summary: ${WEB}"
  done

  echo ""
  log "${PASS} sample(s) with metrics found, ${FAIL} missing."
  log "Multiplex report: ${OUTPUT_DIR}/${RUN_ID}/outs/multiplexing_analysis/"
}

# ─── PRINT SEURAT LOADING CODE ───────────────────────────────────────────────
print_seurat_code() {
  local OUTS_DIR="${OUTPUT_DIR}/${RUN_ID}/outs/per_sample_outs"

  cat <<RCODE

# ==============================================================================
# Load GEM-X Flex 4-plex output into Seurat (R)
# Each sample has its own per_sample_outs/ directory — load and merge.
# ==============================================================================

library(Seurat)
library(dplyr)

outs_root <- "${OUTS_DIR}"
samples    <- c("SC001", "SC004", "SC006", "SC008")

seurat_list <- lapply(samples, function(s) {
  mat_dir <- file.path(outs_root, s, "count", "sample_filtered_feature_bc_matrix")
  counts  <- Read10X(data.dir = mat_dir)
  obj     <- CreateSeuratObject(
               counts    = counts,
               project   = s,
               min.cells = 3,
               min.features = 200
             )
  obj\$sample <- s
  obj
})

# Merge into one object
combined <- merge(
  seurat_list[[1]],
  y        = seurat_list[2:4],
  add.cell.ids = samples,
  project  = "plaque_4plex"
)

# Basic QC metrics
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# QC violin plot per sample
VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  pt.size  = 0,
  ncol     = 3
)

# Typical Flex QC thresholds (adjust to your data)
combined <- subset(
  combined,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt  < 5
)

cat("Cells after QC:", ncol(combined), "\n")
cat("Cells per sample:\n")
print(table(combined\$sample))

RCODE
}

# ─── MAIN ────────────────────────────────────────────────────────────────────
main() {
  log "Pipeline start — GEM-X Flex 4-plex, run: ${RUN_ID}"
  $DRY_RUN && log "[DRY-RUN mode — no commands will be executed]"

  preflight
  CONFIG=$(write_config)
  run_cellranger_multi "$CONFIG"

  if ! $DRY_RUN; then
    qc_summary
  fi

  log "===================================================================="
  log "Output root:        ${OUTPUT_DIR}/${RUN_ID}/outs/"
  log "Per-sample outs:    ${OUTPUT_DIR}/${RUN_ID}/outs/per_sample_outs/"
  log "Multiplex analysis: ${OUTPUT_DIR}/${RUN_ID}/outs/multiplexing_analysis/"
  log "===================================================================="
  log "Seurat loading code:"
  print_seurat_code
  log "Pipeline complete."
}

main
