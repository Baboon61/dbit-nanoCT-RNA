#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

# Keep py_compile cache files out of macOS protected cache locations.
export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-${TMPDIR:-/tmp}/dbit_nanoct_rna_pycache}"

# Fast syntax check for all standalone Python helper scripts.
python3 -m py_compile \
  workflow/scripts/BC_process_CT.py \
  workflow/scripts/BC_process_RNA.py \
  workflow/scripts/check_outputs.py \
  workflow/scripts/debarcode.py \
  workflow/scripts/extract_confident_barcodes.py \
  workflow/scripts/filter_cellranger_gtf_file.py \
  workflow/scripts/remove_LA_duplicates.py \
  workflow/scripts/replace_cellranger_barcodes.py \
  workflow/scripts/restore_cellranger_barcodes.py \
  workflow/scripts/write_pipeline_summary.py \
  workflow/scripts/write_rule_metrics_report.py \
  workflow/scripts/write_tool_version.py

# Build the CT and RNA DAGs without running jobs.
snakemake -n -p -s workflow/Snakefile_CT --use-conda --conda-frontend conda
snakemake -n -p -s workflow/Snakefile_RNA --use-conda --conda-frontend conda
