import argparse
import json
import os
import platform
import subprocess
import sys
from datetime import datetime

try:
    import yaml
except ImportError:
    yaml = None


def command_version(command):
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return None
    first_line = result.stdout.strip().splitlines()
    return first_line[0] if first_line else None


def main():
    parser = argparse.ArgumentParser(description="Write a compact workflow summary and tool-version report.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--workflow", required=True, choices=["CT", "RNA"])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config, "r") as handle:
        config_text = handle.read()

    if yaml is None:
        config = {"raw_config": config_text}
    else:
        config = yaml.safe_load(config_text)

    general = config.get("general", {})
    cellranger = general.get("cellranger", {}) if isinstance(general, dict) else {}

    summary = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "workflow": args.workflow,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "samples": config.get("samples"),
        "general": general,
        "tool_versions": {
            "snakemake": command_version(["snakemake", "--version"]),
            "samtools": command_version(["samtools", "--version"]),
            "bedtools": command_version(["bedtools", "--version"]),
            "macs2": command_version(["macs2", "--version"]),
            "cellranger": command_version([cellranger["software_bin"], "--version"]) if "software_bin" in cellranger else None,
            "fragtk": command_version([general["fragtk_bin"], "--version"]) if "fragtk_bin" in general else None,
        },
    }

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
