# dbit-nanoCT-RNA

Snakemake workflows for processing spatial dbit-style nanoCUT&Tag and RNA sequencing data.

This repository contains two related pipelines:

- `workflow/Snakefile_CT`: spatial nanoCUT&Tag / Cut&Tag processing.
- `workflow/Snakefile_RNA`: spatial RNA processing.

Both workflows share the same early spatial demultiplexing logic, Cell Ranger barcode whitelist handling, path helpers, validation checks, logging, and completion gatekeepers.

## Acknowledgements

This pipeline is based on and adapted from:

- [bartosovic-lab/nanoscope](https://github.com/bartosovic-lab/nanoscope/)
- [dyxmvp/spatial-CUT-Tag](https://github.com/dyxmvp/spatial-CUT-Tag)

The current repository reorganizes those ideas into CT and RNA Snakemake workflows with shared helpers, explicit completion `.done` files, per-rule logging, validation, and additional CT matrix generation.

## Repository Layout

```text
config/
  config_CT.yaml                  # CT sample, barcode, resource, and tool configuration
  config_RNA.yaml                 # RNA sample, resource, and tool configuration

barcodes/
  spatial_barcodes_10000.txt      # Example spatial barcode whitelist
  spatial_barcodes_2500.txt       # Smaller spatial barcode whitelist

envs/
  *.yaml                          # Conda environments used by Snakemake rules

workflow/
  Snakefile_CT                    # Main CT workflow
  Snakefile_RNA                   # Main RNA workflow
  Snakefile_spatial_demux.smk     # Shared primer/linker filtering rules
  Snakefile_cellranger_whitelist.smk
  Snakefile_paths.smk             # Shared path, log, benchmark, and .done helpers
  Snakefile_validation.smk        # Config, input, executable, and directory validation
  Snakefile_func.smk              # FASTQ discovery and filename parsing
  scripts/                        # Python and AWK helper scripts
```

## What The Workflows Do

### Shared CT/RNA Steps

Both workflows:

1. Discover 10X-style FASTQs under `rawData_dir/<sample>/`.
2. Validate config keys, input FASTQ naming, required files, executable paths, and writable directories.
3. Filter reads for the PCR primer, linker 1, and linker 2 with `bbduk.sh`.
4. Build processed FASTQs containing spatial barcode information.
5. Temporarily replace the Cell Ranger barcode whitelist with the configured spatial barcode whitelist.
6. Run Cell Ranger.
7. Restore the original Cell Ranger barcode whitelist.
8. Write logs, benchmarks where configured, `.done` gatekeeper files, and a JSON pipeline summary.

The workflows use `.done` files as completion gatekeepers. A downstream rule depends on the `.done` file from the upstream rule instead of only depending on a partially written data file. This prevents long-running rules from triggering downstream jobs before their outputs are complete.

### CT Workflow

The CT workflow additionally:

1. Splits processed FASTQs by antibody/modality barcode with `debarcode.py`.
2. Runs `cellranger-atac count` separately for each sample/modality/barcode.
3. Removes linear amplification duplicates.
4. Generates no-LA fragment files.
5. Creates normalized bigWig tracks.
6. Calls broad peaks with MACS2.
7. Computes barcode metrics for all reads and peak-overlapping reads.
8. Generates fragment matrices with `fragtk`:
   - MACS broad peak matrix.
   - Fixed genomic bin matrices from `general.bin_sizes`.
   - Gene body plus promoter matrix.

### RNA Workflow

The RNA workflow:

1. Builds the RNA R1 read containing spatial barcode plus UMI.
2. Keeps the sequence read as R2.
3. Runs `cellranger count` using the processed FASTQs.
4. Writes a pipeline summary report.

## Input Requirements

FASTQ files must use 10X-style filenames:

```text
<sample>_<sample_number>_<lane>_<read>_<suffix>.fastq.gz
```

Example:

```text
P35765_1003_S1_L001_R1_001.fastq.gz
P35765_1003_S1_L001_R2_001.fastq.gz
```

Supported extensions:

- `.fastq.gz`
- `.fq.gz`
- `.fastq`
- `.fq`

Expected folder layout:

```text
rawData_dir/
  <sample>/
    ... any nested folders ...
      <sample>_S1_L001_R1_001.fastq.gz
      <sample>_S1_L001_R2_001.fastq.gz
```

`general.read_barcode` controls which read contains the spatial barcode structure:

- `R1`: barcode read is R1 and sequence read is R2.
- `R2`: barcode read is R2 and sequence read is R1.

## Configuration

Edit the relevant config file before running:

- CT: `config/config_CT.yaml`
- RNA: `config/config_RNA.yaml`

### CT Samples

CT samples are configured as a mapping from sample name to antibody/modality barcodes:

```yaml
samples:
  P35765_1003:
    barcodes:
      H3K4me3: CCTATCCT
      H3K27me3: ATAGAGGC
```

Each modality is processed into its own output folder:

```text
processedData_dir/<sample>/<modality>_<barcode>/
```

### RNA Samples

RNA samples are configured as a list:

```yaml
samples:
  - P35765_1010
```

### Important Shared Config Keys

```yaml
general:
  rawData_dir: /path/to/raw
  processedData_dir: /path/to/processed
  spatial_barcodes_file: barcodes/spatial_barcodes_10000.txt
  read_barcode: R2
  core: 8
  tempdir: /path/to/tmp
  cellranger:
    software_bin: /path/to/cellranger-or-cellranger-atac
    reference_path: /path/to/reference/
    barcodes_path: /path/to/cellranger/barcodes/
    barcodes_file: 737K-*.txt.gz
    core: 16
    mem: 64
  PCRprimer_sequence: CAAGCGTTGGCTTCTCGCATCT
  linker1_sequence: GTGGCCGATGTTTCGCATCGGCGTACGACT
  linker2_sequence: ATCCACGTGCTTGAGAGGCCAGAGCATTCG
  spatial_barcode1_length: 8
  spatial_barcode2_length: 8
  hamming_distance_linkers: 3
```

`spatial_barcodes_file` may be absolute or relative to the repository root.

### CT-Specific Config Keys

```yaml
general:
  adapt_Tn5: TAGATGTGTATAAGAGACAG
  nanoCutTag_barcode_length: 8
  fragtk_bin: /path/to/fragtk
  macs_genome: hs
  bin_sizes: [2500]
  debarcode_flush_every: 100000
  debarcode_threads: 1
  debarcode_chunk_size: 10000
  min_first_flush_ratio: 0.00001
```

The CT debarcode rule fails early if the selected antibody barcode is essentially absent in the first flush. This is controlled by:

- `debarcode_flush_every`
- `min_first_flush_ratio`

### RNA-Specific Config Keys

```yaml
general:
  umi_barcode_length: 10
```

## Cell Ranger Barcode Whitelist Handling

Both workflows temporarily replace the Cell Ranger barcode whitelist with the configured spatial barcode list. This allows Cell Ranger to treat spatial barcodes as valid cell barcodes.

The workflow:

1. Backs up the original Cell Ranger barcode whitelist into `processedData_dir/.done/`.
2. Writes a gzipped whitelist generated from `general.spatial_barcodes_file`.
3. Runs Cell Ranger jobs.
4. Restores the original whitelist.

Because this modifies a file inside the Cell Ranger installation, `general.cellranger.barcodes_path` must be writable by the user running Snakemake. Avoid running CT and RNA workflows at the same time if they point to the same Cell Ranger barcode file.

## Installation

Use Conda or Mamba. The base environment contains Snakemake:

```bash
mamba env create -f envs/dbit-spatial-nanoct-rna-base.yaml
conda activate dbit-spatial-nanoct-rna-base
```

Snakemake creates the rule-specific environments from `envs/*.yaml` when run with `--use-conda`.

External tools that must be configured manually:

- Cell Ranger RNA or Cell Ranger ATAC, set in `general.cellranger.software_bin`.
- The matching Cell Ranger reference, set in `general.cellranger.reference_path`.
- `fragtk` for CT matrix generation, set in `general.fragtk_bin`.

## Running

From the repository root:

```bash
snakemake -s workflow/Snakefile_CT --use-conda --cores 32
```

```bash
snakemake -s workflow/Snakefile_RNA --use-conda --cores 32
```

Dry-run:

```bash
snakemake -n -s workflow/Snakefile_CT --use-conda
snakemake -n -s workflow/Snakefile_RNA --use-conda
```

Local workflow check:

```bash
bash workflow/scripts/check_workflow.sh
```

The check script compiles Python helper scripts and runs Snakemake dry-runs for CT and RNA.

## Main Outputs

### Shared Outputs

```text
processedData_dir/
  logs/
  reports/
  .done/
  <sample>/
    logs/
    .done/
    tmp_data/
      qc_raw_data/
      stats/
```

Important shared files:

- `logs/*.log`: per-rule logs.
- `.done/*.done`: completion gatekeepers.
- `reports/pipeline_summary_CT.json` or `reports/pipeline_summary_RNA.json`: compact run metadata and tool versions.
- `all_spatial_cells.txt`: spatial barcode whitelist converted to Cell Ranger-style `<barcode>-1` cell IDs.

### CT Outputs

For each CT sample/modality/barcode:

```text
processedData_dir/<sample>/<modality>_<barcode>/
  fastq/
  cellranger/outs/
    possorted_bam.bam
    possorted_bam.bam.bai
    fragments.tsv.gz
    singlecell.csv
    peaks.bed
    fragments_noLA_duplicates.tsv.gz
    fragments_noLA_duplicates.tsv.gz.tbi
  bigwig/
    all_reads.bw
  peaks/macs_broad/
    <modality>_peaks.broadPeak
  barcode_metrics/
    all_barcodes.txt
    peaks_barcodes.txt
  matrix/
    matrix_peaks/
    matrix_bin_<binsize>/
    matrix_genes/
```

### RNA Outputs

For each RNA sample:

```text
processedData_dir/<sample>/
  cellranger/outs/
    possorted_genome_bam.bam
  tmp_data/
```

## Safeguards And Validation

The workflow fails during DAG creation or immediately after rules when common problems are detected:

- Missing required config keys.
- Invalid integer or ratio parameters.
- Missing raw data directory.
- Missing or empty spatial barcode whitelist.
- Missing or non-executable Cell Ranger / `fragtk` paths.
- Non-writable output, temp, or Cell Ranger barcode directories.
- FASTQ sample prefix mismatch.
- Incomplete R1/R2 groups.
- Wrong FASTQ filename format.
- Empty critical outputs.
- Missing CT barcode matches during the first debarcode flush.

## Implementation Notes

- Path construction is centralized in `workflow/Snakefile_paths.smk`.
- FASTQ discovery is cached during DAG construction in `workflow/Snakefile_func.smk`.
- Shared primer/linker filtering lives in `workflow/Snakefile_spatial_demux.smk`.
- Shared Cell Ranger whitelist replacement/restoration lives in `workflow/Snakefile_cellranger_whitelist.smk`.
- Python scripts are designed to create output directories where needed and to fail clearly when required outputs are missing.
- CT matrix rules write into hidden temporary matrix directories first, validate them, then move them into final locations.

## Development Checks

Run before committing workflow changes:

```bash
bash workflow/scripts/check_workflow.sh
git diff --check
```

If Snakemake is not available in the active shell, at least run:

```bash
PYTHONPYCACHEPREFIX=/tmp/dbit_nanoct_rna_pycache python3 -m py_compile workflow/scripts/*.py
git diff --check
```
