import argparse
import os
import shutil


def main():
    parser = argparse.ArgumentParser(description="Restore the original Cell Ranger barcode whitelist from a backup.")
    parser.add_argument("--cellranger-barcodes", required=True)
    parser.add_argument("--backup", required=True)
    args = parser.parse_args()

    if not os.path.isfile(args.backup):
        raise SystemExit("*** Error: Cell Ranger barcode backup does not exist: {}\n".format(args.backup))

    # Put the Cell Ranger installation back exactly as it was before this run.
    shutil.copy2(args.backup, args.cellranger_barcodes)


if __name__ == "__main__":
    main()
