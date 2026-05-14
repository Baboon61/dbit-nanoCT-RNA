import itertools
import os
import glob
import datetime
import sys
import re

include: 'Snakefile_func.smk'

configfile: workflow.basedir + '/../config/config_RNA.yaml'

#Put CutTag samples into variables
samples_list  = config['samples']
print("RNA samples to be processed : {}".format(samples_list))

#Register if barcodes construction is R1 or R2
if config['general']['read_barcode']=="R1":
  barcode_r1="f"
  barcode_r2="t"
  Rfile_barcode="R1"
  Rfile_sequence="R2"
elif config['general']['read_barcode']=="R2":
  barcode_r1="t"
  barcode_r2="f"
  Rfile_barcode="R2"
  Rfile_sequence="R1"
else:
  print("Error in barcode in R1/R2 in config file, please double check\n")

#Clean path for raw and processed data to have / at the end of the paths
raw_dir = config['general']['rawData_dir'] if config['general']['rawData_dir'][-1] == '/' else config['general']['rawData_dir']+"/"
proc_dir = config['general']['processedData_dir'] if config['general']['processedData_dir'][-1] == '/' else config['general']['processedData_dir']+"/"

def validate_existing_file(path, config_key):
  if not os.path.isfile(path):
    sys.exit("*** Error: {} does not exist or is not a file: {}\n".format(config_key, path))
  if os.path.getsize(path) == 0:
    sys.exit("*** Error: {} is empty: {}\n".format(config_key, path))

def ensure_writable_directory(path, config_key):
  try:
    os.makedirs(path, exist_ok=True)
  except OSError as error:
    sys.exit("*** Error: could not create {} directory '{}': {}\n".format(config_key, path, error))
  if not os.path.isdir(path):
    sys.exit("*** Error: {} is not a directory: {}\n".format(config_key, path))
  if not os.access(path, os.W_OK):
    sys.exit("*** Error: {} is not writable: {}\n".format(config_key, path))

def validate_cellranger_config(cellranger_config):
  if not os.path.isfile(cellranger_config['software_bin']):
    sys.exit("*** Error: general.cellranger.software_bin does not exist or is not a file: {}\n".format(cellranger_config['software_bin']))
  if not os.access(cellranger_config['software_bin'], os.X_OK):
    sys.exit("*** Error: general.cellranger.software_bin is not executable: {}\n".format(cellranger_config['software_bin']))
  if not os.path.isdir(cellranger_config['reference_path']):
    sys.exit("*** Error: general.cellranger.reference_path does not exist or is not a directory: {}\n".format(cellranger_config['reference_path']))
  if not os.path.isdir(cellranger_config['barcodes_path']):
    sys.exit("*** Error: general.cellranger.barcodes_path does not exist or is not a directory: {}\n".format(cellranger_config['barcodes_path']))
  if not os.access(cellranger_config['barcodes_path'], os.W_OK):
    sys.exit("*** Error: general.cellranger.barcodes_path is not writable: {}\n"
             "The workflow replaces the Cell Ranger barcode whitelist in this directory.\n".format(cellranger_config['barcodes_path']))

def validate_sample_fastqs(sample, expected_reads):
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

if not os.path.isdir(raw_dir):
  sys.exit("*** Error: rawData_dir does not exist or is not a directory: {}\n"
           "Please update general.rawData_dir in config/config_RNA.yaml.\n".format(raw_dir))

ensure_writable_directory(proc_dir, "general.processedData_dir")
ensure_writable_directory(config['general']['tempdir'], "general.tempdir")
validate_existing_file(config['general']['spatial_barcodes_file'], "general.spatial_barcodes_file")
validate_cellranger_config(config['general']['cellranger'])

#Check Cellranger format of fastq files
fastq_info_by_sample = {sample_id: validate_sample_fastqs(sample_id, ['R1', 'R2']) for sample_id in samples_list}
res = [res for sample_id in samples_list for res in fastq_info_by_sample[sample_id]]

#Check that the cores argument are lower than the rules cores parameters
if workflow.cores < config['general']['cellranger']['core']:
  sys.exit("*** Error: not enough cores were given to Snakemake.\n"
           "Requested workflow cores: {}\n"
           "Required by general.cellranger.core: {}\n".format(workflow.cores, config['general']['cellranger']['core']))

if workflow.cores < config['general']['core']:
  sys.exit("*** Error: not enough cores were given to Snakemake.\n"
           "Requested workflow cores: {}\n"
           "Required by general.core: {}\n".format(workflow.cores, config['general']['core']))

#Remove the back up file from the cellranger barcodes to allow the copy of the new ones in rule get_barcodes_cellranger
if os.path.exists(config['general']['cellranger']['barcodes_path'] + 'BAK_' + config['general']['cellranger']['barcodes_file']):
  os.remove(config['general']['cellranger']['barcodes_path'] + 'BAK_' + config['general']['cellranger']['barcodes_file'])

#Print the 10X nomenclature error message
Error_message="*** Error: Wrong input files specified. The input must follow 10X fastq files nomenclature, see  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input***\n" +\
        "The files should be placed in the same folder specify as raw folder in the config.yaml file\n"

for sample in samples_list:
  file_name_list = glob.glob(raw_dir + "/" + sample + "/**/*.fastq.gz", recursive=True)
  file_name_list.extend(glob.glob(raw_dir + "/" + sample + "/**/*.fq.gz", recursive=True))
  file_name_list.extend(glob.glob(raw_dir + "/" + sample + "/**/*.fastq", recursive=True))
  file_name_list.extend(glob.glob(raw_dir + "/" + sample + "/**/*.fq", recursive=True))
  for file_name in file_name_list:
    if len(re.split("\.", os.path.basename(file_name))) not in [2, 3]:
      sys.exit(1)
      sys.stderr.write(Error_message)

for file_struc in res:
  if not re.match("S[0-9]+",file_struc['number']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("L[0-9]+",file_struc['lane']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("[RI][0-9]+",file_struc['read']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("[0-9]+",file_struc['suffix']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif file_struc['ext'] not in ['fastq.gz', 'fq.gz', 'fastq', 'fq']:
    sys.exit(1)
    sys.stderr.write(Error_message)
