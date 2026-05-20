import argparse
import os


def main():
    # This helper is called inside Snakemake rules to make empty outputs fatal.
    parser = argparse.ArgumentParser(description="Fail if expected output files are missing or empty.")
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--dirs", nargs="*", default=[])
    args = parser.parse_args()

    # Files must exist and be non-empty; directories must exist and contain files.
    missing = [path for path in args.paths if not os.path.isfile(path)]
    empty = [path for path in args.paths if os.path.isfile(path) and os.path.getsize(path) == 0]
    missing_dirs = [path for path in args.dirs if not os.path.isdir(path)]
    empty_dirs = [path for path in args.dirs if os.path.isdir(path) and not os.listdir(path)]

    if missing or empty or missing_dirs or empty_dirs:
        message = []
        if missing:
            message.append("Missing files:\n  {}".format("\n  ".join(missing)))
        if empty:
            message.append("Empty files:\n  {}".format("\n  ".join(empty)))
        if missing_dirs:
            message.append("Missing directories:\n  {}".format("\n  ".join(missing_dirs)))
        if empty_dirs:
            message.append("Empty directories:\n  {}".format("\n  ".join(empty_dirs)))
        raise SystemExit("*** Error: output validation failed.\n{}\n".format("\n".join(message)))


if __name__ == "__main__":
    main()
