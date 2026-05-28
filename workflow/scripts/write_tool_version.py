import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description="Write the first line of a tool version command.")
    parser.add_argument("--output", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("No version command provided")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        output = result.stdout
    except OSError:
        output = ""
    lines = output.strip().splitlines()
    version = lines[0] if lines else ""

    with open(args.output, "w") as handle:
        handle.write(version)
        handle.write("\n")


if __name__ == "__main__":
    main()
