import argparse
import csv
import gzip
import json
import os
import re
from pathlib import Path

try:
    import yaml
    from yaml.nodes import MappingNode
except ImportError:
    yaml = None
    MappingNode = None

if yaml is not None:
    class MetricsYamlLoader(yaml.SafeLoader):
        pass

    def construct_defaultdict(loader, node):
        if isinstance(node, MappingNode):
            for key_node, value_node in node.value:
                key = loader.construct_object(key_node)
                if key == "dictitems":
                    return dict(loader.construct_mapping(value_node, deep=True))
        return {}

    MetricsYamlLoader.add_constructor(
        "tag:yaml.org,2002:python/object/apply:collections.defaultdict",
        construct_defaultdict,
    )


def read_yaml_config(path):
    with open(path, "r") as handle:
        text = handle.read()
    if yaml is None:
        processed_match = re.search(r"^\s*processedData_dir:\s*(\S+)", text, re.MULTILINE)
        general = {}
        if processed_match:
            general["processedData_dir"] = processed_match.group(1)
        return {"general": general}
    return yaml.safe_load(text) or {}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def count_lines(path):
    try:
        with open_text(path) as handle:
            return sum(1 for _ in handle)
    except OSError:
        return None


def count_bed_records(path):
    try:
        with open_text(path) as handle:
            return sum(1 for line in handle if line.strip() and not line.startswith("#"))
    except OSError:
        return None


def file_size(path):
    try:
        return os.path.getsize(path)
    except OSError:
        return None


def parse_json_file(path):
    try:
        with open(path, "r") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}


def parse_yaml_file(path):
    try:
        with open(path, "r") as handle:
            text = handle.read()
    except OSError:
        return {}
    if yaml is not None:
        data = yaml.load(text, Loader=MetricsYamlLoader)
        return dict(data or {})

    data = {}
    for line in text.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            continue
        try:
            data[key] = int(value)
        except ValueError:
            try:
                data[key] = float(value)
            except ValueError:
                data[key] = value
    return dict(data or {})


def parse_benchmark(path):
    try:
        with open(path, "r") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError:
        return {}
    if not rows:
        return {}
    metrics = {}
    for key, value in rows[-1].items():
        if value is None or value == "":
            continue
        try:
            metrics[key] = float(value)
        except ValueError:
            metrics[key] = value
    return metrics


def parse_bbduk_stats(path):
    metrics = {}
    try:
        with open(path, "r") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = re.split(r"\t+", line)
                if len(fields) >= 2:
                    key = fields[0].strip().lower().replace(" ", "_")
                    try:
                        metrics[key] = int(fields[1])
                    except ValueError:
                        metrics[key] = fields[1]
    except OSError:
        pass
    return metrics


def parse_singlecell_csv(path):
    metrics = {}
    try:
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
    except OSError:
        return metrics

    metrics["barcodes_reported"] = len(rows)
    for column in ["passed_filters", "peak_region_fragments", "TSS_fragments"]:
        if rows and column in rows[0]:
            total = 0
            for row in rows:
                try:
                    total += int(float(row.get(column) or 0))
                except ValueError:
                    pass
            metrics[column + "_sum"] = total

    for column in ["is__cell_barcode", "is_cell_barcode"]:
        if rows and column in rows[0]:
            metrics["cell_barcodes"] = sum(1 for row in rows if str(row.get(column, "")).lower() in {"1", "true", "yes"})
            break
    return metrics


def parse_barcode_counts(path):
    metrics = {"barcode_rows": 0, "read_count_sum": 0}
    try:
        with open(path, "r") as handle:
            for line in handle:
                fields = line.strip().split()
                if len(fields) < 2:
                    continue
                metrics["barcode_rows"] += 1
                try:
                    metrics["read_count_sum"] += int(fields[0])
                except ValueError:
                    pass
    except OSError:
        return {}
    return metrics


def parse_cellranger_metrics_summary(path):
    try:
        with open(path, "r", newline="") as handle:
            rows = list(csv.reader(handle))
    except OSError:
        return {}
    if len(rows) < 2:
        return {}
    metrics = {}
    for key, value in zip(rows[0], rows[1]):
        key = key.strip().lower().replace(" ", "_")
        key = re.sub(r"[^a-z0-9_]+", "", key)
        if key:
            metrics[key] = value.strip()
    return metrics


def matrix_metrics(path):
    files = [item for item in path.rglob("*") if item.is_file()]
    metrics = {
        "file_count": len(files),
        "total_bytes": sum(file_size(item) or 0 for item in files),
    }
    for item in files:
        if item.name.endswith(".mtx") or item.name.endswith(".mtx.gz"):
            metrics["matrix_non_comment_lines"] = count_bed_records(item)
            break
    return metrics


def context_from_benchmark(path, processed_dir):
    rel_parts = path.relative_to(processed_dir).parts
    name = path.stem
    context = {}
    if len(rel_parts) >= 3 and rel_parts[1] == "benchmarks":
        context["sample"] = rel_parts[0]
        context["scope"] = "sample"
    elif len(rel_parts) >= 4 and rel_parts[2] == "benchmarks":
        context["sample"] = rel_parts[0]
        modality_barcode = rel_parts[1]
        if "_" in modality_barcode:
            context["modality"], context["barcode"] = modality_barcode.rsplit("_", 1)
        context["scope"] = "modality"
    elif len(rel_parts) >= 2 and rel_parts[0] == "benchmarks":
        context["scope"] = "global"
    else:
        context["scope"] = "unknown"

    rule = name
    for prefix in [
        "filter_primer", "filter_L1", "filter_L2", "bc_process", "seq_file_rename",
        "debarcode", "matrix_bin",
    ]:
        if name.startswith(prefix):
            rule = prefix
            break
    rule = {
        "macs_broad": "run_macs_broad",
        "noLA_fragments": "bam_noLA_to_fragments_noLA",
        "matrix_peaks": "create_matrix_peaks",
        "matrix_genes": "create_genebody_and_promoter_matrix",
        "matrix_bin": "create_matrix_bins",
    }.get(rule, rule)
    context["rule"] = rule
    context["benchmark_file"] = str(path)
    return context


def target_label(rule, workflow):
    labels = {
        "CT": {
            "filter_primer": "filterPCR",
            "filter_L1": "filterL1",
            "filter_L2": "filterL2",
        },
        "RNA": {
            "filter_primer": "filterprimer",
            "filter_L1": "filterL1",
            "filter_L2": "filterL2",
        },
    }
    return labels.get(workflow, {}).get(rule, rule)


def add_entry(entries, rule, context, workflow, metrics=None, benchmark=None, outputs=None):
    entry = {
        "rule": rule,
        "target_label": target_label(rule, workflow),
        "context": context,
        "runtime": benchmark or {},
        "metrics": metrics or {},
        "outputs": outputs or [],
    }
    entries.append(entry)


def collect_benchmark_entries(processed_dir, workflow):
    entries = {}
    for path in processed_dir.rglob("benchmarks/*.txt"):
        context = context_from_benchmark(path, processed_dir)
        rule = context.pop("rule")
        key = (rule, context.get("sample"), context.get("modality"), context.get("barcode"), path.stem)
        entries[key] = {
            "rule": rule,
            "target_label": target_label(rule, workflow),
            "context": context,
            "runtime": parse_benchmark(path),
            "metrics": {},
            "outputs": [],
        }
    return entries


def attach_or_add(entries, rule, context, workflow, metrics, outputs):
    for entry in entries.values():
        if entry["rule"] != rule:
            continue
        entry_context = entry["context"]
        if all(entry_context.get(key) == value for key, value in context.items()):
            entry["metrics"].update(metrics)
            entry["outputs"].extend(outputs)
            return
    key = (rule, context.get("sample"), context.get("modality"), context.get("barcode"), len(entries))
    entries[key] = {
        "rule": rule,
        "target_label": target_label(rule, workflow),
        "context": context,
        "runtime": {},
        "metrics": metrics,
        "outputs": outputs,
    }


def collect_common_output_metrics(processed_dir, workflow):
    entries = collect_benchmark_entries(processed_dir, workflow)

    for path in processed_dir.glob("*/tmp_data/stats/*_bc_process_*.json"):
        sample = path.parts[-4]
        attach_or_add(entries, "bc_process", {"sample": sample}, workflow, parse_json_file(path), [str(path)])

    for path in processed_dir.glob("*/tmp_data/qc_raw_data/*_stats.*.txt"):
        sample = path.parts[-4]
        step_match = re.search(r"_stats\.(.+)\.txt$", path.name)
        step = step_match.group(1) if step_match else "unknown"
        rule = {"primer": "filter_primer", "linker1": "filter_L1", "linker2": "filter_L2"}.get(step, step)
        attach_or_add(entries, rule, {"sample": sample}, workflow, parse_bbduk_stats(path), [str(path)])

    return entries


def collect_ct_output_metrics(processed_dir, workflow):
    entries = collect_common_output_metrics(processed_dir, workflow)

    for path in processed_dir.glob("*/*_*/fastq/*_statistics.yaml"):
        sample = path.parts[-4]
        modality_barcode = path.parts[-3]
        context = {"sample": sample}
        if "_" in modality_barcode:
            context["modality"], context["barcode"] = modality_barcode.rsplit("_", 1)
        stats = parse_yaml_file(path)
        metrics = {"debarcoded_reads": sum(value for value in stats.values() if isinstance(value, int))}
        metrics.update(stats)
        attach_or_add(entries, "debarcode", context, workflow, metrics, [str(path)])

    for modality_dir in processed_dir.glob("*/*_*"):
        if not modality_dir.is_dir():
            continue
        sample = modality_dir.parts[-2]
        modality_barcode = modality_dir.name
        context = {"sample": sample}
        if "_" in modality_barcode:
            context["modality"], context["barcode"] = modality_barcode.rsplit("_", 1)

        outs = modality_dir / "cellranger" / "outs"
        if outs.is_dir():
            metrics = {}
            singlecell = outs / "singlecell.csv"
            if singlecell.exists():
                metrics.update(parse_singlecell_csv(singlecell))
            fragments = outs / "fragments.tsv.gz"
            if fragments.exists():
                metrics["fragments"] = count_lines(fragments)
            peaks = outs / "peaks.bed"
            if peaks.exists():
                metrics["cellranger_peaks"] = count_bed_records(peaks)
            bam = outs / "possorted_bam.bam"
            if bam.exists():
                metrics["possorted_bam_bytes"] = file_size(bam)
            attach_or_add(entries, "run_cellranger", context, workflow, metrics, [str(item) for item in [singlecell, fragments, peaks, bam] if item.exists()])

        for stats_file in (modality_dir / "cellranger" / "outs").glob("*_stats.txt"):
            attach_or_add(entries, "remove_LA_duplicates", context, workflow, parse_json_file(stats_file), [str(stats_file)])

        no_la_fragments = modality_dir / "cellranger" / "outs" / "fragments_noLA_duplicates.tsv.gz"
        if no_la_fragments.exists():
            attach_or_add(entries, "sort_sinto_output", context, workflow, {"noLA_fragments": count_lines(no_la_fragments)}, [str(no_la_fragments)])

        broad_peak = modality_dir / "peaks" / "macs_broad"
        for peak_file in broad_peak.glob("*_peaks.broadPeak"):
            attach_or_add(entries, "run_macs_broad", context, workflow, {"macs_broad_peaks": count_bed_records(peak_file)}, [str(peak_file)])

        bigwig_dir = modality_dir / "bigwig"
        for bigwig in bigwig_dir.glob("*.bw"):
            attach_or_add(entries, "bam_to_bw", context, workflow, {"bigwig_bytes": file_size(bigwig)}, [str(bigwig)])

        barcode_dir = modality_dir / "barcode_metrics"
        for counts in barcode_dir.glob("*.txt"):
            rule = "barcode_metrics_peaks" if counts.name.startswith("peaks_") else "barcode_metrics_all"
            attach_or_add(entries, rule, context, workflow, parse_barcode_counts(counts), [str(counts)])

        matrix_dir = modality_dir / "matrix"
        if matrix_dir.is_dir():
            for child in matrix_dir.iterdir():
                if child.is_dir() and not child.name.startswith("."):
                    if child.name.startswith("matrix_bin"):
                        rule = "create_matrix_bins"
                    elif child.name == "matrix_peaks":
                        rule = "create_matrix_peaks"
                    elif child.name == "matrix_genes":
                        rule = "create_genebody_and_promoter_matrix"
                    else:
                        rule = "create_" + child.name
                    attach_or_add(entries, rule, context | {"matrix": child.name}, workflow, matrix_metrics(child), [str(child)])

    return sorted(entries.values(), key=lambda item: (item["rule"], json.dumps(item["context"], sort_keys=True)))


def collect_rna_output_metrics(processed_dir, workflow):
    entries = collect_common_output_metrics(processed_dir, workflow)

    for sample_dir in processed_dir.iterdir() if processed_dir.exists() else []:
        if not sample_dir.is_dir() or sample_dir.name.startswith("."):
            continue
        context = {"sample": sample_dir.name}
        outs = sample_dir / "cellranger" / "outs"
        if not outs.is_dir():
            continue

        metrics = {}
        outputs = []
        metrics_summary = outs / "metrics_summary.csv"
        if metrics_summary.exists():
            metrics.update(parse_cellranger_metrics_summary(metrics_summary))
            outputs.append(str(metrics_summary))

        bam = outs / "possorted_genome_bam.bam"
        if bam.exists():
            metrics["possorted_genome_bam_bytes"] = file_size(bam)
            outputs.append(str(bam))

        for matrix_name in ["filtered_feature_bc_matrix", "raw_feature_bc_matrix"]:
            matrix_dir = outs / matrix_name
            if not matrix_dir.is_dir():
                continue
            matrix = matrix_metrics(matrix_dir)
            barcodes = matrix_dir / "barcodes.tsv.gz"
            features = matrix_dir / "features.tsv.gz"
            matrix_file = matrix_dir / "matrix.mtx.gz"
            if barcodes.exists():
                matrix[matrix_name + "_barcodes"] = count_lines(barcodes)
            if features.exists():
                matrix[matrix_name + "_features"] = count_lines(features)
            if matrix_file.exists():
                matrix[matrix_name + "_matrix_lines"] = count_lines(matrix_file)
            metrics.update(matrix)
            outputs.append(str(matrix_dir))

        molecule_info = outs / "molecule_info.h5"
        if molecule_info.exists():
            metrics["molecule_info_h5_bytes"] = file_size(molecule_info)
            outputs.append(str(molecule_info))

        attach_or_add(entries, "run_cellranger", context, workflow, metrics, outputs)

    return sorted(entries.values(), key=lambda item: (item["rule"], json.dumps(item["context"], sort_keys=True)))


def main():
    parser = argparse.ArgumentParser(description="Write a rule runtime and output metrics report.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--workflow", required=True, choices=["CT", "RNA"])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml_config(args.config)
    processed_dir = Path(config.get("general", {}).get("processedData_dir", "."))
    if args.workflow == "CT":
        entries = collect_ct_output_metrics(processed_dir, args.workflow)
    else:
        entries = collect_rna_output_metrics(processed_dir, args.workflow)

    report = {
        "workflow": args.workflow,
        "processedData_dir": str(processed_dir),
        "rules": entries,
    }

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
