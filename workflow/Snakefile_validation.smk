import os
import re
import sys

# Resolve paths written relative to the repository root, for example barcodes/*.txt.
def resolve_repo_path(path):
  return path if os.path.isabs(path) else os.path.join(repo_root, path)

# Fail during DAG creation if a required file is missing or empty.
def validate_existing_file(path, config_key):
  if not os.path.isfile(path):
    sys.exit("*** Error: {} does not exist or is not a file: {}\n".format(config_key, path))
  if os.path.getsize(path) == 0:
    sys.exit("*** Error: {} is empty: {}\n".format(config_key, path))

# Executables such as Cell Ranger and fragtk must exist and be runnable.
def validate_executable(path, config_key):
  validate_existing_file(path, config_key)
  if not os.access(path, os.X_OK):
    sys.exit("*** Error: {} is not executable: {}\n".format(config_key, path))

def require_config_keys(container, required_keys, config_key):
  missing_keys = [key for key in required_keys if key not in container]
  if missing_keys:
    sys.exit("*** Error: missing required config key(s) under {}: {}\n".format(
      config_key,
      ", ".join(missing_keys)
    ))

def validate_positive_int(value, config_key):
  if not isinstance(value, int) or value <= 0:
    sys.exit("*** Error: {} must be a positive integer, got: {}\n".format(config_key, value))

def validate_nonnegative_int(value, config_key):
  if not isinstance(value, int) or value < 0:
    sys.exit("*** Error: {} must be a non-negative integer, got: {}\n".format(config_key, value))

def validate_ratio(value, config_key):
  if not isinstance(value, (int, float)) or value < 0 or value > 1:
    sys.exit("*** Error: {} must be a number between 0 and 1, got: {}\n".format(config_key, value))

def with_trailing_slash(path):
  return path if path[-1] == '/' else path + '/'

def get_read_layout():
  if config['general']['read_barcode'] == "R1":
    return "f", "t", "R1", "R2"
  if config['general']['read_barcode'] == "R2":
    return "t", "f", "R2", "R1"
  sys.exit("*** Error: general.read_barcode must be either R1 or R2.\n")

def validate_general_config_schema(modality):
  required_general = [
    'rawData_dir',
    'processedData_dir',
    'spatial_barcodes_file',
    'read_barcode',
    'core',
    'tempdir',
    'cellranger',
    'PCRprimer_sequence',
    'linker1_sequence',
    'linker2_sequence',
    'spatial_barcode1_length',
    'spatial_barcode2_length',
    'hamming_distance_linkers'
  ]
  require_config_keys(config, ['samples'], 'root')
  require_config_keys(config['general'], required_general, 'general')
  require_config_keys(config['general']['cellranger'], ['software_bin', 'reference_path', 'barcodes_path', 'barcodes_file', 'core', 'mem'], 'general.cellranger')

  if config['general']['read_barcode'] not in ['R1', 'R2']:
    sys.exit("*** Error: general.read_barcode must be either R1 or R2.\n")

  validate_positive_int(config['general']['core'], 'general.core')
  validate_positive_int(config['general']['cellranger']['core'], 'general.cellranger.core')
  validate_positive_int(config['general']['cellranger']['mem'], 'general.cellranger.mem')
  validate_positive_int(config['general']['spatial_barcode1_length'], 'general.spatial_barcode1_length')
  validate_positive_int(config['general']['spatial_barcode2_length'], 'general.spatial_barcode2_length')
  validate_nonnegative_int(config['general']['hamming_distance_linkers'], 'general.hamming_distance_linkers')

  if modality == 'CT':
    require_config_keys(config['general'], ['fragtk_bin', 'adapt_Tn5', 'nanoCutTag_barcode_length', 'macs_genome', 'bin_sizes'], 'general')
    validate_positive_int(config['general']['nanoCutTag_barcode_length'], 'general.nanoCutTag_barcode_length')
    for key in ['samtools_threads', 'sinto_threads', 'bamcoverage_threads', 'debarcode_threads', 'debarcode_chunk_size', 'debarcode_flush_every']:
      if key in config['general']:
        validate_positive_int(config['general'][key], 'general.{}'.format(key))
    for key in ['samtools_mem_mb', 'sinto_mem_mb', 'bamcoverage_mem_mb', 'macs_mem_mb', 'remove_LA_duplicates_mem_mb']:
      if key in config['general']:
        validate_positive_int(config['general'][key], 'general.{}'.format(key))
    if 'min_first_flush_ratio' in config['general']:
      validate_ratio(config['general']['min_first_flush_ratio'], 'general.min_first_flush_ratio')
    if not isinstance(config['general']['bin_sizes'], list) or not config['general']['bin_sizes']:
      sys.exit("*** Error: general.bin_sizes must be a non-empty list of positive integers.\n")
    for binsize in config['general']['bin_sizes']:
      validate_positive_int(binsize, 'general.bin_sizes')

  if modality == 'RNA':
    require_config_keys(config['general'], ['umi_barcode_length'], 'general')
    validate_positive_int(config['general']['umi_barcode_length'], 'general.umi_barcode_length')

# Output and temp directories are created early so failures happen before jobs run.
def ensure_writable_directory(path, config_key):
  try:
    os.makedirs(path, exist_ok=True)
  except OSError as error:
    sys.exit("*** Error: could not create {} directory '{}': {}\n".format(config_key, path, error))
  if not os.path.isdir(path):
    sys.exit("*** Error: {} is not a directory: {}\n".format(config_key, path))
  if not os.access(path, os.W_OK):
    sys.exit("*** Error: {} is not writable: {}\n".format(config_key, path))

# Cell Ranger config is special because the workflow temporarily replaces its whitelist.
def validate_cellranger_config(cellranger_config):
  validate_executable(cellranger_config['software_bin'], "general.cellranger.software_bin")
  if not os.path.isdir(cellranger_config['reference_path']):
    sys.exit("*** Error: general.cellranger.reference_path does not exist or is not a directory: {}\n".format(cellranger_config['reference_path']))
  if not os.path.isdir(cellranger_config['barcodes_path']):
    sys.exit("*** Error: general.cellranger.barcodes_path does not exist or is not a directory: {}\n".format(cellranger_config['barcodes_path']))
  if not os.access(cellranger_config['barcodes_path'], os.W_OK):
    sys.exit("*** Error: general.cellranger.barcodes_path is not writable: {}\n"
             "The workflow needs temporary write access to replace and restore the Cell Ranger barcode whitelist.\n".format(cellranger_config['barcodes_path']))
  validate_existing_file(
    cellranger_config['barcodes_path'] + cellranger_config['barcodes_file'],
    "general.cellranger.barcodes_path/general.cellranger.barcodes_file"
  )

# Validate that each sample has complete 10X-style FASTQ groups.
def validate_sample_fastqs(sample, expected_reads, raw_dir, config_path):
  records = get_fastq_info_from_folder(raw_dir, sample)
  if not records:
    sys.exit("*** Error: no FASTQ files were found for sample: {}\n"
             "Expected files under rawData_dir/<sample>/ matching **/*.fastq.gz, **/*.fq.gz, **/*.fastq, or **/*.fq.\n"
             "rawData_dir is currently: {}\n".format(sample, raw_dir))

  unexpected_ids = sorted(set(record['id'] for record in records if record['id'] != sample))
  if unexpected_ids:
    sys.exit("*** Error: FASTQ sample ID mismatch for configured sample '{}'.\n"
             "Found file prefix(es): {}\n"
             "FASTQ names must start with the configured sample name, for example '{}_S1_L001_R1_001.fastq.gz'.\n"
             .format(sample, ", ".join(unexpected_ids), sample))

  read_groups = {}
  for record in records:
    key = (record['id'], record['number'], record['lane'], record['suffix'], record['ext'])
    read_groups.setdefault(key, set()).add(record['read'])

  missing_groups = []
  for key, reads in sorted(read_groups.items()):
    missing_reads = sorted(set(expected_reads) - reads)
    if missing_reads:
      missing_groups.append("{} {} {} missing {}".format(key[1], key[2], key[3], ",".join(missing_reads)))

  if missing_groups:
    sys.exit("*** Error: incomplete FASTQ read pairs for sample '{}'.\n"
             "Each sample/lane/suffix group must contain {}.\n"
             "Problems:\n  {}\n".format(sample, ", ".join(expected_reads), "\n  ".join(missing_groups)))

  return records

# parse_fastq() already checks the filename shape; this keeps explicit field checks.
def validate_fastq_records(records):
  for record in records:
    if not re.match("S[0-9]+", record['number']):
      sys.exit("*** Error: wrong FASTQ sample number: {}\n".format(record['number']))
    if not re.match("L[0-9]+", record['lane']):
      sys.exit("*** Error: wrong FASTQ lane: {}\n".format(record['lane']))
    if not re.match("[RI][0-9]+", record['read']):
      sys.exit("*** Error: wrong FASTQ read: {}\n".format(record['read']))
    if not re.match("[0-9]+", record['suffix']):
      sys.exit("*** Error: wrong FASTQ suffix: {}\n".format(record['suffix']))
    if record['ext'] not in ['fastq.gz', 'fq.gz', 'fastq', 'fq']:
      sys.exit("*** Error: wrong FASTQ extension: {}\n".format(record['ext']))

# Common CT/RNA config validation used before Snakemake builds the full DAG.
def validate_common_config(config_path):
  if not os.path.isdir(raw_dir):
    sys.exit("*** Error: rawData_dir does not exist or is not a directory: {}\n"
             "Please update general.rawData_dir in {}.\n".format(raw_dir, config_path))

  ensure_writable_directory(proc_dir, "general.processedData_dir")
  ensure_writable_directory(config['general']['tempdir'], "general.tempdir")
  validate_existing_file(config['general']['spatial_barcodes_file'], "general.spatial_barcodes_file")
  validate_cellranger_config(config['general']['cellranger'])
