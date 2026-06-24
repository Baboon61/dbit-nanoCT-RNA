import argparse
import csv
import html
import json
import os
from collections import defaultdict


DEFAULT_RULE_ORDER = [
    "filter_primer",
    "filter_L1",
    "filter_L2",
    "bc_process",
    "debarcode",
    "get_barcodes_cellranger",
    "run_cellranger",
    "remove_LA_duplicates",
    "sort_sinto_output",
    "bam_to_bw",
    "run_macs_broad",
    "barcode_metrics_peaks",
    "barcode_metrics_all",
    "create_matrix_peaks",
    "create_matrix_bins",
    "create_genebody_and_promoter_matrix",
]


RULE_ALIASES = {
    "macs_broad": "run_macs_broad",
    "matrix_peaks": "create_matrix_peaks",
    "matrix_genes": "create_genebody_and_promoter_matrix",
    "remove_LA_duplicates": "remove_LA_duplicates",
    "possort_noLA_bam": "possort_noLA_bam_file",
    "noLA_fragments": "bam_noLA_to_fragments_noLA",
    "sort_noLA_fragments": "sort_sinto_output",
    "clean_cellranger": "clean_cellranger_output",
}


HIDDEN_CONTEXT_KEYS = {"benchmark_file", "scope"}
HIDDEN_RULES = {
    "all_spatial_cells",
    "bam_noLA_to_fragments_noLA",
    "bam_to_namesorted",
    "clean_cellranger_output",
    "pipeline_summary",
    "possort_noLA_bam_file",
    "rule_metrics_html",
    "rule_metrics_report",
    "seq_file_rename",
    "spatial_barcodes_to_cells",
}
RUNTIME_KEYS = ["h:m:s"]
CONTEXT_ORDER = ["sample", "modality", "barcode", "number", "lane", "suffix", "ext", "matrix"]
GOOD_KEYWORDS = ("kept", "matched", "passed", "written", "cell_barcodes", "fragments", "peaks")
BAD_KEYWORDS = ("discarded", "failed", "duplicates")
FILTER_RULES = {"filter_primer", "filter_L1", "filter_L2"}
FILTER_HIDDEN_METRICS = {
    "1_percent",
    "1_reads",
    "input_reads",
    "raw_reads_start",
    "reads_discarded",
    "reads_discarded_percent_vs_previous",
}
BC_PROCESS_HIDDEN_METRICS = {
    "reads_discarded",
    "reads_discarded_percent_vs_previous",
    "reads_processed",
}
DEBARCODE_HIDDEN_METRICS = {
    "debarcoded_reads",
    "no_barcode_match",
    "no_spacer_found",
    "too_short_read",
}
RUN_CELLRANGER_HIDDEN_METRICS = {"cell_barcodes", "fragment_count_histogram", "fragments"}
RUN_CELLRANGER_FIRST_METRICS = [
    "passed_filters_sum",
    "cellranger_peaks",
    "peak_region_fragments_sum",
    "TSS_fragments_sum",
    "median_tss_enrichment_score",
]
BARCODE_METRICS_RULES = {"barcode_metrics_peaks", "barcode_metrics_all"}
BARCODE_METRICS_HIDDEN_METRICS = {"barcode_rows"}
MATRIX_RULES = {"create_matrix_peaks", "create_matrix_bins", "create_genebody_and_promoter_matrix"}
RULE_METRIC_LABELS = {
    ("barcode_metrics_all", "read_count_sum"): "total read count",
    ("barcode_metrics_peaks", "median_read_count"): "median read count in peaks",
    ("create_genebody_and_promoter_matrix", "features_tsv_lines"): "number of genes",
    ("create_matrix_bins", "features_tsv_lines"): "number of bins",
    ("create_matrix_peaks", "features_tsv_lines"): "number of peaks",
}
RULE_METRIC_DESCRIPTIONS = {
    ("barcode_metrics_all", "median_read_count"): "Median of the read-count column in barcode_metrics/all_barcodes.txt.",
    ("barcode_metrics_all", "read_count_sum"): "Sum of the read-count column in barcode_metrics/all_barcodes.txt.",
    ("barcode_metrics_peaks", "median_read_count"): "Median of the read-count column in barcode_metrics/peaks_barcodes.txt for reads overlapping peaks.",
    ("barcode_metrics_peaks", "read_count_sum"): "Sum of the read-count column in barcode_metrics/peaks_barcodes.txt for reads overlapping peaks.",
    ("create_genebody_and_promoter_matrix", "features_tsv_lines"): "Number of feature rows in matrix/matrix_genes/features.tsv.gz.",
    ("create_matrix_bins", "features_tsv_lines"): "Number of feature rows in each matrix/matrix_bin_<binsize>/features.tsv.gz file.",
    ("create_matrix_peaks", "features_tsv_lines"): "Number of feature rows in matrix/matrix_peaks/features.tsv.gz.",
}
METRIC_LABELS = {
    "passed_filters_sum": "pass QC fragments",
    "peak_region_fragments_sum": "fragments in peaks",
    "bigwig_bytes": "bigwig file size",
    "features_tsv_lines": "features.tsv.gz lines",
    "median_read_count": "median read count",
    "possorted_bam_bytes": "bam file size",
    "read_count_sum": "total read count in peaks",
    "TSS_fragments_sum": "fragments in TSS",
    "median_tss_enrichment_score": "median TSS enrichment score",
    "tss_enrichment_score": "median TSS enrichment score",
}
RULE_LABELS = {
    "bc_process": "BC process",
}
METRIC_DESCRIPTIONS = {
    "1_percent": "Percentage value reported by the upstream filtering statistics file.",
    "1_reads": "Read count reported by the upstream filtering statistics file.",
    "barcodes_reported": "Number of rows in Cell Ranger singlecell.csv after excluding the NO_BARCODE row.",
    "bigwig_bytes": "File size of the generated bigWig file.",
    "cell_barcodes": "Number of Cell Ranger singlecell.csv rows flagged as cell barcodes.",
    "cellranger_peaks": "Number of non-empty, non-comment records in Cell Ranger outs/peaks.bed.",
    "file_count": "Number of files found in this output directory.",
    "fragments": "Line count of Cell Ranger outs/fragments.tsv.gz.",
    "fragment_count_histogram": "Histogram data used to draw the fragments-per-cell graph.",
    "input_reads": "Reads entering this rule or filtering step.",
    "LA duplicates": "Read pairs with the same read 1 position and cell barcode but a different mate position; these are removed from the no-LA BAM.",
    "macs_broad_peaks": "Number of non-empty, non-comment records in the MACS broadPeak output.",
    "matrix_non_comment_lines": "Number of non-comment records in the matrix file.",
    "median_read_count": "Median of the read-count column in the barcode metrics file.",
    "median_tss_enrichment_score": "Median TSS enrichment score reported by Cell Ranger summary.csv or metrics_summary.csv.",
    "molecule_info_h5_bytes": "File size of Cell Ranger molecule_info.h5.",
    "noLA_fragments": "Sum of the last column in fragments_noLA_duplicates.tsv.gz after LA duplicate removal.",
    "no_match": "Debarcoding reads not assigned to the selected barcode after no-barcode, no-spacer, and too-short categories are collapsed for display.",
    "not mapped": "Read pairs skipped during LA duplicate removal because at least one mate is unmapped.",
    "not proper": "Read pairs skipped during LA duplicate removal because at least one mate is not marked as a proper pair.",
    "passed_filters_sum": "Sum of Cell Ranger singlecell.csv passed_filters across reported barcode rows, excluding NO_BARCODE.",
    "PCR duplicates": "Read pairs with the same read 1 position, mate position, and cell barcode as a previously observed pair; these are retained in the no-LA BAM.",
    "peak_region_fragments_sum": "Sum of Cell Ranger singlecell.csv peak_region_fragments across reported barcode rows, excluding NO_BARCODE.",
    "possorted_bam_bytes": "File size of Cell Ranger outs/possorted_bam.bam.",
    "possorted_genome_bam_bytes": "File size of Cell Ranger outs/possorted_genome_bam.bam.",
    "raw_reads_start": "Raw reads entering the first primer-filtering step.",
    "read_count_sum": "Sum of the read-count column in the barcode metrics file.",
    "reads_discarded": "Reads removed by the current processing step.",
    "reads_discarded_percent_vs_previous": "Discarded reads divided by input reads for this step.",
    "reads_kept": "Reads written by the current filtering or processing step.",
    "reads_kept_percent_vs_previous": "Reads kept divided by input reads for this step.",
    "reads_kept_percent_vs_raw": "Reads kept divided by raw reads at the start of the workflow.",
    "reads_processed": "Reads processed by the barcode-processing script.",
    "reads_selected": "Reads assigned to the selected antibody/modality barcode during debarcoding.",
    "reads_selected_percent_vs_bc_process": "Selected reads divided by the reads written by BC process for the same sample.",
    "reads_written": "Reads written by the barcode-processing script.",
    "total_bytes": "Total file size across files in this output directory.",
    "TSS_fragments_sum": "Sum of Cell Ranger singlecell.csv TSS_fragments across reported barcode rows, excluding NO_BARCODE.",
    "tss_enrichment": "TSS enrichment value reported by Cell Ranger.",
    "tss_enrichment_score": "TSS enrichment score reported by Cell Ranger; displayed as median TSS enrichment score when no median-specific field exists.",
    "unique": "First observed read-pair position for a cell barcode during LA duplicate removal.",
}


def escape(value):
    return html.escape("" if value is None else str(value))


def human_label(value):
    if value in METRIC_LABELS:
        return escape(METRIC_LABELS[value])
    if value in RULE_LABELS:
        return escape(RULE_LABELS[value])
    words = ["TSS" if word.lower() == "tss" else word for word in str(value).replace("_", " ").split(" ")]
    return escape(" ".join(words))


def metric_label(key, rule=None):
    if (rule, key) in RULE_METRIC_LABELS:
        return escape(RULE_METRIC_LABELS[(rule, key)])
    return human_label(key)


def metric_description(key, rule=None):
    if (rule, key) in RULE_METRIC_DESCRIPTIONS:
        return RULE_METRIC_DESCRIPTIONS[(rule, key)]
    if key in METRIC_DESCRIPTIONS:
        return METRIC_DESCRIPTIONS[key]
    return f"Metric generated by {human_label(rule) if rule else 'this rule'} from key '{key}'."


def render_metric_header(key, rule=None):
    description = escape(metric_description(key, rule=rule))
    return (
        '<span class="metric-label">'
        f'<span>{metric_label(key, rule=rule)}</span>'
        f'<span class="metric-help" tabindex="0" title="{description}" data-tooltip="{description}" aria-label="{description}">!</span>'
        '</span>'
    )


def css_suffix(value):
    return "".join(char if char.isalnum() or char in {"-", "_"} else "-" for char in str(value))


def format_value(value):
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if abs(value) >= 100:
            return f"{value:,.1f}"
        if abs(value) >= 1:
            return f"{value:,.3f}"
        return f"{value:.3g}"
    if isinstance(value, (list, tuple)):
        return ", ".join(format_value(item) for item in value)
    if isinstance(value, dict):
        return escape(json.dumps(value, sort_keys=True))
    return escape(value)


def format_bytes(value):
    try:
        size = float(value)
    except (TypeError, ValueError):
        return format_value(value)
    units = ["B", "KB", "MB", "GB", "TB"]
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:,.1f} {unit}" if unit != "B" else f"{int(size):,} {unit}"
        size /= 1024
    return format_value(value)


def numeric_value(value):
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip().rstrip("%"))
        except ValueError:
            return None
    return None


def format_metric(key, value):
    if key.endswith("_bytes") or key == "bytes":
        return format_bytes(value)
    if "percent" in key.lower() or key.lower().endswith("pct"):
        value_number = numeric_value(value)
        if value_number is not None:
            return f"{value_number:,.2f}%"
    return format_value(value)


def metric_status(key, value):
    key_lower = key.lower()
    value_number = numeric_value(value)

    if key_lower == "raw_reads_start" or key_lower == "input_reads":
        return None
    if "percent" in key_lower or key_lower.endswith("pct"):
        if value_number is None:
            return None
        if any(word in key_lower for word in BAD_KEYWORDS):
            if value_number <= 20:
                return "good"
            if value_number <= 50:
                return "warn"
            return "bad"
        if any(word in key_lower for word in GOOD_KEYWORDS):
            if value_number >= 80:
                return "good"
            if value_number >= 50:
                return "warn"
            return "bad"
        return None
    if any(word in key_lower for word in BAD_KEYWORDS):
        return "bad"
    if any(word in key_lower for word in GOOD_KEYWORDS):
        return "good"
    return None


def retention_color(value):
    value_number = numeric_value(value)
    if value_number is None:
        return None
    red = (249, 216, 214)
    yellow = (255, 241, 199)
    green = (220, 239, 235)
    if value_number < 70:
        color = red
    elif value_number > 95:
        color = green
    elif value_number <= 82.5:
        ratio = (value_number - 70) / 12.5
        color = tuple(round(red[index] + (yellow[index] - red[index]) * ratio) for index in range(3))
    else:
        ratio = (value_number - 82.5) / 12.5
        color = tuple(round(yellow[index] + (green[index] - yellow[index]) * ratio) for index in range(3))
    return "#{:02x}{:02x}{:02x}".format(*color)


def debarcode_selection_color(value):
    value_number = numeric_value(value)
    if value_number is None:
        return None
    red = (249, 216, 214)
    yellow = (255, 241, 199)
    green = (220, 239, 235)
    def blend(start, end, ratio):
        return tuple(round(start[index] + (end[index] - start[index]) * ratio) for index in range(3))
    if 45 <= value_number <= 55:
        color = green
    elif value_number < 35 or value_number > 65:
        color = red
    elif value_number < 45:
        if value_number <= 40:
            color = blend(red, yellow, (value_number - 35) / 5)
        else:
            color = blend(yellow, green, (value_number - 40) / 5)
    else:
        if value_number <= 60:
            color = blend(green, yellow, (value_number - 55) / 5)
        else:
            color = blend(yellow, red, (value_number - 60) / 5)
    return "#{:02x}{:02x}{:02x}".format(*color)


def render_metric_value(key, value, values=None, rule=None):
    formatted = format_metric(key, value)
    if rule == "remove_LA_duplicates" and key in {"not mapped", "not proper"}:
        return f'<span class="metric-tag warn">{formatted}</span>'
    if rule == "remove_LA_duplicates" and key == "unique":
        return f'<span class="metric-tag good">{formatted}</span>'
    if rule in BARCODE_METRICS_RULES and key in {"read_count_sum", "median_read_count"}:
        return f'<span class="metric-tag good">{formatted}</span>'
    if rule in MATRIX_RULES and key == "features_tsv_lines":
        return f'<span class="metric-tag good">{formatted}</span>'
    if key in {"barcodes_reported", "bigwig_bytes", "possorted_bam_bytes"}:
        return f'<span class="metric-tag info">{formatted}</span>'
    if key == "median_tss_enrichment_score":
        return f'<span class="metric-tag good">{formatted}</span>'
    if key == "no_match":
        return f'<span class="metric-tag warn">{formatted}</span>'
    if key == "reads_kept_percent_vs_previous":
        color = retention_color(value)
        if color:
            return f'<span class="metric-tag" style="background-color: {color}; color: #18202a;">{formatted}</span>'
    if key == "reads_selected_percent_vs_bc_process":
        color = debarcode_selection_color(value)
        if color:
            return f'<span class="metric-tag" style="background-color: {color}; color: #18202a;">{formatted}</span>'
    if key in {"reads_kept_percent_vs_raw", "reads_kept", "reads_written"} and values:
        color = retention_color(values.get("reads_kept_percent_vs_previous"))
        if color:
            return f'<span class="metric-tag" style="background-color: {color}; color: #18202a;">{formatted}</span>'
    if key == "reads_selected" and values:
        color = debarcode_selection_color(values.get("reads_selected_percent_vs_bc_process"))
        if color:
            return f'<span class="metric-tag" style="background-color: {color}; color: #18202a;">{formatted}</span>'
    status = metric_status(key, value)
    if not status:
        return formatted
    return f'<span class="metric-tag {status}">{formatted}</span>'


def normalize_rule(rule):
    if rule.startswith("matrix_bin_"):
        return "create_matrix_bins"
    return RULE_ALIASES.get(rule, rule)


def clean_context(entry):
    rule = normalize_rule(entry.get("rule", "unknown"))
    context = {
        key: value
        for key, value in (entry.get("context") or {}).items()
        if key not in HIDDEN_CONTEXT_KEYS
    }

    benchmark_file = (entry.get("context") or {}).get("benchmark_file", "")
    benchmark_name = os.path.splitext(os.path.basename(benchmark_file))[0]
    if rule == "create_matrix_bins" and benchmark_name.startswith("matrix_bin"):
        context.setdefault("matrix", benchmark_name)
    elif rule == "create_matrix_peaks" and benchmark_name == "matrix_peaks":
        context.setdefault("matrix", "matrix_peaks")
    elif rule == "create_genebody_and_promoter_matrix" and benchmark_name == "matrix_genes":
        context.setdefault("matrix", "matrix_genes")

    return context


def entry_context_label(entry):
    context = clean_context(entry)
    parts = []
    for key in CONTEXT_ORDER:
        if key in context:
            parts.append(f"{key}: {context[key]}")
    for key in sorted(set(context) - set(CONTEXT_ORDER)):
        parts.append(f"{key}: {context[key]}")
    return "; ".join(parts) or "global"


def render_context_tags(entry):
    context = clean_context(entry)
    if not context:
        return ""
    keys = [key for key in CONTEXT_ORDER if key in context]
    keys.extend(key for key in sorted(set(context) - set(CONTEXT_ORDER)))
    tags = []
    for key in keys:
        tags.append(
            f'<span class="context-tag context-tag-{css_suffix(key)}">'
            f'<span>{human_label(key)}</span>'
            f'<b>{escape(context[key])}</b>'
            '</span>'
        )
    return '<div class="context-tags">' + "".join(tags) + "</div>"


def context_key(context):
    return tuple((key, context[key]) for key in CONTEXT_ORDER if key in context) + tuple(
        (key, context[key]) for key in sorted(set(context) - set(CONTEXT_ORDER))
    )


def contexts_match(left, right):
    common_keys = set(left) & set(right)
    if not common_keys:
        return not left and not right
    return all(left[key] == right[key] for key in common_keys)


def merge_outputs(left, right):
    merged = list(left or [])
    seen = set(merged)
    for item in right or []:
        if item not in seen:
            merged.append(item)
            seen.add(item)
    return merged


def sanitize_runtime(runtime):
    return {key: runtime[key] for key in RUNTIME_KEYS if key in (runtime or {})}


def merge_entries(entries):
    merged = []
    for entry in entries:
        normalized = {
            "rule": normalize_rule(entry.get("rule", "unknown")),
            "target_label": entry.get("target_label") or normalize_rule(entry.get("rule", "unknown")),
            "context": clean_context(entry),
            "runtime": sanitize_runtime(entry.get("runtime") or {}),
            "metrics": dict(entry.get("metrics") or {}),
            "outputs": list(entry.get("outputs") or []),
        }

        match = None
        for candidate in merged:
            if candidate["rule"] == normalized["rule"] and contexts_match(candidate["context"], normalized["context"]):
                match = candidate
                break

        if match is None:
            merged.append(normalized)
            continue

        match["context"].update(normalized["context"])
        match["metrics"].update(normalized["metrics"])
        match["runtime"].update(normalized["runtime"])
        match["outputs"] = merge_outputs(match["outputs"], normalized["outputs"])

    return sorted(merged, key=lambda item: (item["rule"], context_key(item["context"])))


def render_key_values(values, class_name="kv", hidden_keys=None, first_keys=None, last_keys=None, rule=None):
    hidden_keys = hidden_keys or set()
    first_keys = first_keys or []
    last_keys = last_keys or set()
    first_keys = [key for key in first_keys if key in values and key not in hidden_keys]
    visible_keys = [
        key for key in sorted(values)
        if key not in hidden_keys and key not in first_keys and key not in last_keys
    ]
    visible_keys = first_keys + visible_keys
    visible_keys.extend(key for key in sorted(last_keys) if key in values and key not in hidden_keys)
    if not visible_keys:
        return '<span class="muted">none</span>'
    rows = []
    for key in visible_keys:
        rows.append(
            "<tr>"
            f"<th>{render_metric_header(key, rule=rule)}</th>"
            f"<td>{render_metric_value(key, values[key], values, rule=rule)}</td>"
            "</tr>"
        )
    return f'<table class="{class_name}"><tbody>' + "".join(rows) + "</tbody></table>"


def relative_output_path(path, processed_dir):
    if not processed_dir or not os.path.isabs(str(path)):
        return path
    try:
        return os.path.relpath(path, processed_dir)
    except ValueError:
        return path


def render_outputs(outputs, processed_dir=None):
    if not outputs:
        return '<span class="muted">none</span>'
    items = "".join(f"<li>{escape(relative_output_path(path, processed_dir))}</li>" for path in outputs)
    return f"<ul>{items}</ul>"


def output_path(path, processed_dir=None):
    if os.path.isabs(str(path)) or not processed_dir:
        return str(path)
    return os.path.join(processed_dir, str(path))


def metric_summary_key(value):
    normalized = str(value).strip().lower().replace(" ", "_")
    return "".join(char for char in normalized if char.isalnum() or char == "_")


def read_tss_enrichment(path):
    try:
        with open(path, "r", newline="") as handle:
            rows = list(csv.reader(handle))
    except OSError:
        return None
    if len(rows) < 2:
        return None
    summary = {metric_summary_key(key): value.strip() for key, value in zip(rows[0], rows[1])}
    for key in ["median_tss_enrichment_score", "median_tss_enrichment", "tss_enrichment_score", "tss_enrichment"]:
        if summary.get(key):
            return summary[key]
    return None


def add_tss_enrichment_metric(entry, processed_dir=None):
    if entry.get("rule") != "run_cellranger":
        return
    metrics = entry.setdefault("metrics", {})
    if metrics.get("median_tss_enrichment_score") is not None:
        metrics.pop("tss_enrichment_score", None)
        metrics.pop("tss_enrichment", None)
        return
    candidates = []
    for path in entry.get("outputs") or []:
        full_path = output_path(path, processed_dir)
        basename = os.path.basename(full_path)
        if basename in {"metrics_summary.csv", "summary.csv"}:
            candidates.append(full_path)
        elif basename == "singlecell.csv":
            candidates.append(os.path.join(os.path.dirname(full_path), "metrics_summary.csv"))
            candidates.append(os.path.join(os.path.dirname(full_path), "summary.csv"))
    for path in candidates:
        value = read_tss_enrichment(path)
        if value is not None:
            metrics["median_tss_enrichment_score"] = value
            metrics.pop("tss_enrichment_score", None)
            metrics.pop("tss_enrichment", None)
            return
    for key in ["tss_enrichment_score", "tss_enrichment"]:
        if metrics.get(key) is not None:
            metrics["median_tss_enrichment_score"] = metrics.pop(key)
            return


def histogram_cell_count(item):
    if not isinstance(item, dict):
        return None
    value = numeric_value(item.get("cells"))
    return int(value) if value is not None else None


def histogram_x_value(item, key):
    if not isinstance(item, dict):
        return None
    value = numeric_value(item.get(key))
    if value is not None:
        return value
    label = str(item.get("label", ""))
    if "-" in label:
        label = label.split("-", 1)[0 if key == "start" else 1]
    return numeric_value(label)


def trim_histogram_plateau(histogram):
    items = [item for item in histogram if isinstance(item, dict) and histogram_cell_count(item) is not None]
    if len(items) < 8:
        return items
    counts = [histogram_cell_count(item) for item in items]
    max_count = max(counts)
    tail_count = counts[-1]
    tail_start = len(counts) - 1
    while tail_start > 0 and counts[tail_start - 1] == tail_count:
        tail_start -= 1
    tail_length = len(counts) - tail_start
    if tail_length >= 4 and tail_count <= max(2, max_count * 0.1):
        return items[: tail_start + 1]
    return items


def rounded_secondary_tick(value):
    if value is None:
        return None
    value = int(value)
    if value <= 0:
        return 0
    if value >= 1000:
        unit = 1000
    elif value >= 100:
        unit = 100
    elif value >= 10:
        unit = 10
    else:
        unit = 1
    return max(unit, (value // unit) * unit)


def axis_ticks(maximum, count=5):
    if maximum is None:
        return []
    maximum = int(round(maximum))
    if count <= 1 or maximum <= 0:
        return [0]
    ticks = [0]
    for index in range(1, count - 1):
        ticks.append(rounded_secondary_tick(maximum * index / (count - 1)))
    ticks.append(maximum)
    return ticks


def format_axis_tick(value):
    if value is None:
        return ""
    return format_value(int(round(value)))


def render_fragment_histogram(entry):
    if entry.get("rule") != "run_cellranger":
        return ""
    histogram = (entry.get("metrics") or {}).get("fragment_count_histogram") or []
    if not isinstance(histogram, list):
        return ""
    histogram = trim_histogram_plateau(histogram)
    bars = []
    cell_counts = [histogram_cell_count(item) for item in histogram]
    cell_counts = [count for count in cell_counts if count is not None]
    max_cells = max(cell_counts, default=0)
    if max_cells <= 0:
        return ""
    y_ticks = []
    for value in reversed(axis_ticks(max_cells, count=3)):
        bottom = 0 if max_cells == 0 else max(0, min(100, (value / max_cells) * 100))
        y_ticks.append(
            f'<div class="histogram-y-tick" style="bottom: {bottom}%;">'
            f'<span>{format_axis_tick(value)}</span>'
            '</div>'
        )

    for item in histogram:
        cells = histogram_cell_count(item) or 0
        label = str(item.get("label", ""))
        height = max(3, round((cells / max_cells) * 100))
        title = f'{label} fragments: {format_value(int(cells))} cells'
        bars.append(
            '<div class="histogram-bar" '
            f'style="height: {height}%;" title="{escape(title)}"></div>'
        )
    if not bars:
        return ""
    max_x = histogram_x_value(histogram[-1], "end")
    x_ticks = "".join(f"<span>{escape(format_axis_tick(value))}</span>" for value in axis_ticks(max_x, count=5))
    return (
        '<section class="fragment-graph">'
        '<h4>Fragments per cell</h4>'
        '<div class="histogram-area">'
        '<div class="histogram-y-ticks">'
        + "".join(y_ticks)
        + "</div>"
        '<div class="histogram" aria-label="Cell count by fragment count">'
        + "".join(bars)
        + "</div>"
        + "</div>"
        + f'<div class="histogram-x-ticks">{x_ticks}</div>'
        "</section>"
    )


def render_entry(entry, processed_dir=None):
    add_tss_enrichment_metric(entry, processed_dir=processed_dir)
    target = entry.get("target_label") or entry.get("rule")
    context = render_context_tags(entry)
    runtime = entry.get("runtime") or {}
    runtime_tag = f'<span class="time-tag">{escape(runtime["h:m:s"])}</span>' if runtime.get("h:m:s") else ""
    hidden_metrics = FILTER_HIDDEN_METRICS if entry.get("rule") in FILTER_RULES else set()
    if entry.get("rule") == "bc_process":
        hidden_metrics = BC_PROCESS_HIDDEN_METRICS
    if entry.get("rule") == "debarcode":
        hidden_metrics = DEBARCODE_HIDDEN_METRICS
    if entry.get("rule") == "run_cellranger":
        hidden_metrics = RUN_CELLRANGER_HIDDEN_METRICS
    if entry.get("rule") in BARCODE_METRICS_RULES:
        hidden_metrics = BARCODE_METRICS_HIDDEN_METRICS
    last_metrics = {"no_match"} if entry.get("rule") == "debarcode" else set()
    first_metrics = []
    if entry.get("rule") == "run_cellranger":
        first_metrics = RUN_CELLRANGER_FIRST_METRICS
        last_metrics = {"barcodes_reported", "possorted_bam_bytes"}
    metrics = render_key_values(
        entry.get("metrics") or {},
        hidden_keys=hidden_metrics,
        first_keys=first_metrics,
        last_keys=last_metrics,
        rule=entry.get("rule"),
    )
    graph = render_fragment_histogram(entry)
    outputs = render_outputs(entry.get("outputs") or [], processed_dir=processed_dir)
    return f"""
      <article class="entry">
        <div class="entry-head">
          <div class="entry-title">
            <h3>{escape(target)}</h3>
            {runtime_tag}
          </div>
          {context}
        </div>
        <div class="entry-grid">
          <section>
            <h4>Metrics</h4>
            {metrics}
          </section>
          {graph}
          <section class="outputs">
            <h4>Outputs</h4>
            {outputs}
          </section>
        </div>
      </article>
    """


def rule_sort_key(rule, order_index):
    normalized = normalize_rule(rule)
    return (order_index.get(normalized, len(order_index)), normalized)


def total_raw_reads(entries):
    total = 0
    found = False
    for entry in entries:
        if entry.get("rule") != "filter_primer":
            continue
        metrics = entry.get("metrics") or {}
        value = metrics.get("raw_reads_start", metrics.get("input_reads"))
        if value is None:
            continue
        try:
            total += int(value)
            found = True
        except (TypeError, ValueError):
            continue
    return total if found else None


def build_html(report, rule_order):
    all_entries = merge_entries(report.get("rules") or [])
    entries = [entry for entry in all_entries if entry.get("rule") not in HIDDEN_RULES]
    order_index = {rule: index for index, rule in enumerate(rule_order)}
    groups = defaultdict(list)
    for entry in entries:
        groups[entry.get("rule", "unknown")].append(entry)

    raw_reads = total_raw_reads(all_entries)
    raw_reads_label = f"{raw_reads:,}" if raw_reads is not None else "n/a"
    workflow = report.get("workflow", "workflow")
    processed_dir = report.get("processedData_dir", "")

    sections = []
    for rule in sorted(groups, key=lambda item: rule_sort_key(item, order_index)):
        rule_entries = sorted(groups[rule], key=lambda item: entry_context_label(item))
        sections.append(
            f"""
            <section class="rule-section">
              <div class="rule-title">
                <h2>{human_label(rule)}</h2>
                <span>{len(rule_entries):,} entr{"y" if len(rule_entries) == 1 else "ies"}</span>
              </div>
              {''.join(render_entry(entry, processed_dir=processed_dir) for entry in rule_entries)}
            </section>
            """
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(workflow)} rule metrics</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f7f9;
      --panel: #ffffff;
      --text: #18202a;
      --muted: #657080;
      --line: #dbe1e8;
      --accent: #136f63;
      --accent-soft: #dcefeb;
      --shadow: 0 1px 2px rgba(24, 32, 42, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.45;
    }}
    header {{
      background: var(--panel);
      border-bottom: 1px solid var(--line);
      padding: 28px clamp(18px, 4vw, 48px);
    }}
    h1, h2, h3, h4, p {{ margin: 0; }}
    h1 {{
      font-size: clamp(1.8rem, 3vw, 2.6rem);
      font-weight: 700;
    }}
    .subtitle {{
      color: var(--muted);
      margin-top: 8px;
      max-width: 900px;
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 24px clamp(14px, 3vw, 32px) 48px;
    }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin-bottom: 24px;
    }}
    .stat {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px 16px;
      box-shadow: var(--shadow);
    }}
    .stat b {{
      display: block;
      font-size: 1.55rem;
      line-height: 1.1;
    }}
    .stat span {{
      color: var(--muted);
      font-size: 0.9rem;
    }}
    .rule-section {{
      margin-top: 20px;
      border-top: 3px solid var(--accent);
      padding-top: 12px;
    }}
    .rule-title {{
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 16px;
      margin-bottom: 10px;
    }}
    .rule-title h2 {{
      font-size: 1.25rem;
      text-transform: capitalize;
    }}
    .rule-title span {{
      background: #e8edf3;
      border-radius: 999px;
      color: #354052;
      padding: 3px 10px;
      font-size: 0.86rem;
      font-weight: 650;
      white-space: nowrap;
    }}
    .entry {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      box-shadow: var(--shadow);
      margin: 10px 0;
      overflow: hidden;
    }}
    .entry-head {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
      background: #fbfcfd;
    }}
    .entry-title {{
      align-items: center;
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }}
    .entry-head h3 {{
      font-size: 1rem;
      font-weight: 650;
    }}
    .time-tag, .context-tag {{
      border-radius: 999px;
      display: inline-flex;
      line-height: 1.2;
      white-space: nowrap;
    }}
    .time-tag {{
      background: #e8edf3;
      color: #354052;
      font-size: 0.82rem;
      font-variant-numeric: tabular-nums;
      font-weight: 650;
      padding: 3px 9px;
    }}
    .context-tags {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      justify-content: flex-end;
    }}
    .context-tag {{
      background: #eef7f6;
      border: 1px solid #c7e0dc;
      color: #184e48;
      font-size: 0.9rem;
      overflow: hidden;
    }}
    .context-tag span {{
      background: #dcefeb;
      color: #42645f;
      padding: 3px 7px;
    }}
    .context-tag b {{
      font-weight: 650;
      padding: 3px 8px;
    }}
    .context-tag-sample {{
      background: #eaf2ff;
      border-color: #bfd4f5;
      color: #173b73;
    }}
    .context-tag-sample span {{
      background: #d7e7ff;
      color: #315f9d;
    }}
    .context-tag-modality {{
      background: #f2ecff;
      border-color: #d7c6f4;
      color: #4b2c78;
    }}
    .context-tag-modality span {{
      background: #e5d8fb;
      color: #664996;
    }}
    .context-tag-barcode {{
      background: #fff3df;
      border-color: #f2c985;
      color: #714600;
    }}
    .context-tag-barcode span {{
      background: #ffe4b8;
      color: #8a5c12;
    }}
    .context-tag-lane {{
      background: #e8f6ee;
      border-color: #b9dcc8;
      color: #195737;
    }}
    .context-tag-lane span {{
      background: #d5efdf;
      color: #31704b;
    }}
    .entry-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: 0;
    }}
    .entry-grid section {{
      padding: 12px 14px;
      border-right: 1px solid var(--line);
      min-width: 0;
    }}
    .entry-grid section:last-child {{
      border-right: 0;
    }}
    h4 {{
      color: var(--muted);
      font-size: 0.78rem;
      letter-spacing: 0;
      margin-bottom: 8px;
      text-transform: uppercase;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.9rem;
    }}
    th, td {{
      border-top: 1px solid #edf1f4;
      padding: 5px 0;
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-weight: 500;
      padding-right: 12px;
      text-align: left;
      width: 52%;
    }}
    .metric-label {{
      align-items: center;
      display: inline-flex;
      gap: 5px;
    }}
    .metric-help {{
      align-items: center;
      background: #e8edf3;
      border-radius: 50%;
      color: #354052;
      cursor: help;
      display: inline-flex;
      font-size: 0.62rem;
      font-weight: 800;
      height: 14px;
      justify-content: center;
      line-height: 1;
      position: relative;
      width: 14px;
    }}
    .metric-help:hover::after,
    .metric-help:focus::after {{
      background: #18202a;
      border-radius: 4px;
      bottom: calc(100% + 6px);
      color: #ffffff;
      content: attr(data-tooltip);
      font-size: 0.72rem;
      font-weight: 500;
      left: 50%;
      line-height: 1.3;
      max-width: min(320px, 70vw);
      min-width: 220px;
      padding: 7px 8px;
      position: absolute;
      transform: translateX(-50%);
      white-space: normal;
      z-index: 5;
    }}
    .metric-help:focus {{
      outline: 2px solid #9fc0ee;
      outline-offset: 2px;
    }}
    td {{
      font-variant-numeric: tabular-nums;
      overflow-wrap: anywhere;
    }}
    ul {{
      margin: 0;
      padding-left: 18px;
      color: var(--muted);
      font-size: 0.86rem;
      overflow-wrap: anywhere;
    }}
    .fragment-graph {{
      display: flex;
      flex-direction: column;
      gap: 6px;
    }}
    .histogram-area {{
      display: grid;
      gap: 6px;
      grid-template-columns: 34px minmax(0, 1fr);
      margin-top: auto;
    }}
    .histogram-y-ticks {{
      color: var(--muted);
      font-size: 0.68rem;
      position: relative;
    }}
    .histogram-y-tick {{
      left: 0;
      position: absolute;
      right: 0;
    }}
    .histogram-y-tick span {{
      display: block;
      line-height: 1;
      transform: translateY(-50%);
    }}
    .histogram {{
      align-items: flex-end;
      border-bottom: 1px solid var(--line);
      border-left: 1px solid var(--line);
      display: flex;
      gap: 2px;
      height: 120px;
      min-width: 0;
      overflow-x: auto;
      padding: 8px 4px 0 8px;
    }}
    .histogram-bar {{
      background: #d7e7ff;
      border-radius: 3px 3px 0 0;
      flex: 1 0 8px;
      min-width: 8px;
    }}
    .histogram-x-ticks {{
      color: var(--muted);
      display: flex;
      font-size: 0.68rem;
      justify-content: space-between;
      margin-bottom: auto;
      padding-left: 40px;
    }}
    .metric-tag {{
      display: inline-block;
      border-radius: 999px;
      font-weight: 650;
      line-height: 1.2;
      padding: 2px 8px;
    }}
    .metric-tag.good {{
      background: #dcefeb;
      color: #0d5f55;
    }}
    .metric-tag.warn {{
      background: #fff1c7;
      color: #7a4e00;
    }}
    .metric-tag.bad {{
      background: #f9d8d6;
      color: #9b2c24;
    }}
    .metric-tag.info {{
      background: #d7e7ff;
      color: #173b73;
    }}
    .muted {{
      color: var(--muted);
      font-size: 0.9rem;
    }}
    @media (max-width: 860px) {{
      .entry-head {{
        display: block;
      }}
      .context-tags {{
        justify-content: flex-start;
        margin-top: 8px;
      }}
      .entry-grid {{
        grid-template-columns: 1fr;
      }}
      .entry-grid section {{
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }}
      .entry-grid section:last-child {{
        border-bottom: 0;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>{escape(workflow)} rule metrics</h1>
    <p class="subtitle">Metrics from : {escape(processed_dir)}</p>
  </header>
  <main>
    <section class="summary">
      <div class="stat"><b>{len(all_entries):,}</b><span>rule entries</span></div>
      <div class="stat"><b>{raw_reads_label}</b><span>raw reads at start</span></div>
    </section>
    {''.join(sections)}
  </main>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Render rule metrics JSON as a readable HTML report.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rule-order", nargs="*", default=DEFAULT_RULE_ORDER)
    args = parser.parse_args()

    with open(args.input, "r") as handle:
        report = json.load(handle)

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write(build_html(report, args.rule_order))


if __name__ == "__main__":
    main()
