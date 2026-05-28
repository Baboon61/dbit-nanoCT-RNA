import argparse
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

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


def recorded_tool_versions(paths):
    versions = {}
    for path in paths:
        tool = Path(path).parent.name
        try:
            with open(path, "r") as handle:
                version = handle.readline().strip()
        except OSError:
            continue
        if version and tool not in versions:
            versions[tool] = version
    return versions


def main():
    parser = argparse.ArgumentParser(description="Write a compact workflow summary and tool-version report.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--workflow", required=True, choices=["CT", "RNA"])
    parser.add_argument("--output", required=True)
    parser.add_argument("--tool-version-files", nargs="*", default=[])
    args = parser.parse_args()

    with open(args.config, "r") as handle:
        config_text = handle.read()

    if yaml is None:
        config = {"raw_config": config_text}
    else:
        config = yaml.safe_load(config_text)

    general = config.get("general", {})
    cellranger = general.get("cellranger", {}) if isinstance(general, dict) else {}
    tool_versions = recorded_tool_versions(args.tool_version_files)

    summary = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "workflow": args.workflow,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "samples": config.get("samples"),
        "general": general,
        "tool_versions": {
            "snakemake": tool_versions.get("snakemake") or command_version(["snakemake", "--version"]),
            "samtools": tool_versions.get("samtools") or command_version(["samtools", "--version"]),
            "bedtools": tool_versions.get("bedtools") or command_version(["bedtools", "--version"]),
            "bamcoverage": tool_versions.get("bamcoverage") or command_version(["bamCoverage", "--version"]),
            "macs2": tool_versions.get("macs2") or command_version(["macs2", "--version"]),
            "sinto": tool_versions.get("sinto") or command_version(["sinto", "--version"]),
            "cellranger": tool_versions.get("cellranger") or (command_version([cellranger["software_bin"], "--version"]) if "software_bin" in cellranger else None),
            "fragtk": tool_versions.get("fragtk") or (command_version([general["fragtk_bin"], "--version"]) if "fragtk_bin" in general else None),
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
