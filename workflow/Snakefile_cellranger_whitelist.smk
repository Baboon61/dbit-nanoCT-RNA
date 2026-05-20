# Shared rules that temporarily replace Cell Ranger's barcode whitelist.
# CT and RNA both need Cell Ranger to see the spatial barcodes as valid cells.

rule get_barcodes_cellranger:
  input:
    spatial_barcodes = config['general']['spatial_barcodes_file'],
  output:
    backup = cellranger_barcodes_backup(),
    outfile = touch(global_done('get_barcodes_cellranger'))
  log:
    global_log('get_barcodes_cellranger')
  params:
    script = workflow.basedir + '/scripts/replace_cellranger_barcodes.py',
    cellranger_barcodes_file = cellranger_barcodes_path()
  shell:
    '''
    python3 {params.script} \
      --spatial-barcodes {input.spatial_barcodes} \
      --cellranger-barcodes {params.cellranger_barcodes_file} \
      --backup {output.backup} > {log} 2>&1
    '''

rule restore_cellranger_barcodes:
  input:
    backup = cellranger_barcodes_backup(),
    run_cellranger = cellranger_done_targets
  output:
    outfile = touch(global_done('restore_cellranger_barcodes'))
  log:
    global_log('restore_cellranger_barcodes')
  params:
    script = workflow.basedir + '/scripts/restore_cellranger_barcodes.py',
    cellranger_barcodes_file = cellranger_barcodes_path()
  shell:
    '''
    python3 {params.script} \
      --cellranger-barcodes {params.cellranger_barcodes_file} \
      --backup {input.backup} > {log} 2>&1
    '''
