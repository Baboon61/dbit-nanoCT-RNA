import argparse
import gzip
import os
import shutil
import tempfile


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def write_gzip_from_lines(lines, output_path):
    # Write atomically: create a temp gzip in the target directory, then rename it.
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix=".barcodes.", suffix=".tmp.gz", dir=output_dir)
    os.close(fd)
    try:
        with gzip.open(tmp_path, "wt") as out:
            for line in lines:
                barcode = line.strip().split()[0]
                if barcode:
                    out.write(barcode + "\n")
        os.replace(tmp_path, output_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def main():
    parser = argparse.ArgumentParser(description="Replace a Cell Ranger barcode whitelist with a run-specific spatial whitelist.")
    parser.add_argument("--spatial-barcodes", required=True)
    parser.add_argument("--cellranger-barcodes", required=True)
    parser.add_argument("--backup", required=True)
    args = parser.parse_args()

    backup_dir = os.path.dirname(args.backup)
    if backup_dir:
        os.makedirs(backup_dir, exist_ok=True)

    if not os.path.exists(args.backup):
        # Keep one workflow-local copy of the original whitelist for restoration.
        shutil.copy2(args.cellranger_barcodes, args.backup)
    else:
        # If rerunning after a previous replacement, restore the original first.
        shutil.copy2(args.backup, args.cellranger_barcodes)

    # Cell Ranger expects a gzipped one-column whitelist file.
    with open_text(args.spatial_barcodes) as spatial_barcodes:
        write_gzip_from_lines(spatial_barcodes, args.cellranger_barcodes)


if __name__ == "__main__":
    main()
