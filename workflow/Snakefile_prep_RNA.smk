import os
import sys

include: 'Snakefile_func.smk'

configfile: workflow.basedir + '/../config/config_RNA.yaml'

repo_root = os.path.abspath(os.path.join(workflow.basedir, '..'))

include: 'Snakefile_validation.smk'

#Put RNA samples into variables
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
  sys.exit("*** Error: general.read_barcode must be either R1 or R2.\n")

#Clean path for raw and processed data to have / at the end of the paths
raw_dir = config['general']['rawData_dir'] if config['general']['rawData_dir'][-1] == '/' else config['general']['rawData_dir']+"/"
proc_dir = config['general']['processedData_dir'] if config['general']['processedData_dir'][-1] == '/' else config['general']['processedData_dir']+"/"
config['general']['spatial_barcodes_file'] = resolve_repo_path(config['general']['spatial_barcodes_file'])

validate_common_config("config/config_RNA.yaml")

#Check Cellranger format of fastq files
fastq_info_by_sample = {sample_id: validate_sample_fastqs(sample_id, ['R1', 'R2'], raw_dir, "config/config_RNA.yaml") for sample_id in samples_list}
res = [res for sample_id in samples_list for res in fastq_info_by_sample[sample_id]]
validate_fastq_records(res)
