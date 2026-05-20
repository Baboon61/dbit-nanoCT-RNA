import itertools
import os
import sys

include: 'Snakefile_func.smk'

configfile: workflow.basedir + '/../config/config_CT.yaml'

repo_root = os.path.abspath(os.path.join(workflow.basedir, '..'))

include: 'Snakefile_validation.smk'

#Put CutTag samples into variables
samples_list  = list(config['samples'].keys())
print("CutTag samples to be processed : {}".format(samples_list))

#Put barcodes into variables
barcodes_dict = {sample: config['samples'][sample]['barcodes'] for sample in samples_list}
print("List of barcodes to be processed : {}".format(barcodes_dict))

#Put antibodies into variables
antibodies_list = list(set(itertools.chain(*[barcodes_dict[sample].keys() for sample in barcodes_dict.keys()])))
print("List of antibodies to be processed : {}".format(antibodies_list))

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

validate_common_config("config/config_CT.yaml")
validate_executable(config['general']['fragtk_bin'], "general.fragtk_bin")

#Check Cellranger format of fastq files
fastq_info_by_sample = {sample_id: validate_sample_fastqs(sample_id, ['R1', 'R2'], raw_dir, "config/config_CT.yaml") for sample_id in samples_list}
res = [res for sample_id in samples_list for res in fastq_info_by_sample[sample_id]]
validate_fastq_records(res)
