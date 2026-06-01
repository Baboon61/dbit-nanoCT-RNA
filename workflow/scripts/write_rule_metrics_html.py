import argparse
import html
import json
import os
from collections import defaultdict


DEFAULT_RULE_ORDER = [
    "filter_primer",
    "filter_L1",
    "filter_L2",
    "bc_process",
    "seq_file_rename",
    "debarcode",
    "get_barcodes_cellranger",
    "run_cellranger",
    "bam_to_namesorted",
    "remove_LA_duplicates",
    "possort_noLA_bam_file",
    "bam_noLA_to_fragments_noLA",
    "sort_sinto_output",
    "bam_to_bw",
    "run_macs_broad",
    "barcode_metrics_peaks",
    "barcode_metrics_all",
    "clean_cellranger_output",
    "spatial_barcodes_to_cells",
    "create_matrix_peaks",
    "create_matrix_bins",
    "create_genebody_and_promoter_matrix",
    "pipeline_summary",
    "rule_metrics_report",
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


def escape(value):
    return html.escape("" if value is None else str(value))


def human_label(value):
    return escape(str(value).replace("_", " "))


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


def render_metric_value(key, value, values=None):
    formatted = format_metric(key, value)
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
        return '<span class="muted">global</span>'
    keys = [key for key in CONTEXT_ORDER if key in context]
    keys.extend(key for key in sorted(set(context) - set(CONTEXT_ORDER)))
    tags = []
    for key in keys:
        tags.append(
            '<span class="context-tag">'
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


def render_key_values(values, class_name="kv", hidden_keys=None, last_keys=None):
    hidden_keys = hidden_keys or set()
    last_keys = last_keys or set()
    visible_keys = [key for key in sorted(values) if key not in hidden_keys and key not in last_keys]
    visible_keys.extend(key for key in sorted(last_keys) if key in values and key not in hidden_keys)
    if not visible_keys:
        return '<span class="muted">none</span>'
    rows = []
    for key in visible_keys:
        rows.append(
            "<tr>"
            f"<th>{human_label(key)}</th>"
            f"<td>{render_metric_value(key, values[key], values)}</td>"
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


def render_entry(entry, processed_dir=None):
    target = entry.get("target_label") or entry.get("rule")
    context = render_context_tags(entry)
    runtime = entry.get("runtime") or {}
    runtime_tag = f'<span class="time-tag">{escape(runtime["h:m:s"])}</span>' if runtime.get("h:m:s") else ""
    hidden_metrics = FILTER_HIDDEN_METRICS if entry.get("rule") in FILTER_RULES else set()
    if entry.get("rule") == "bc_process":
        hidden_metrics = BC_PROCESS_HIDDEN_METRICS
    if entry.get("rule") == "debarcode":
        hidden_metrics = DEBARCODE_HIDDEN_METRICS
    last_metrics = {"no_match"} if entry.get("rule") == "debarcode" else set()
    metrics = render_key_values(entry.get("metrics") or {}, hidden_keys=hidden_metrics, last_keys=last_metrics)
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
    entries = merge_entries(report.get("rules") or [])
    order_index = {rule: index for index, rule in enumerate(rule_order)}
    groups = defaultdict(list)
    for entry in entries:
        groups[entry.get("rule", "unknown")].append(entry)

    raw_reads = total_raw_reads(entries)
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
      color: var(--accent);
      background: var(--accent-soft);
      border-radius: 999px;
      padding: 3px 10px;
      font-size: 0.86rem;
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
    .entry-grid {{
      display: grid;
      grid-template-columns: minmax(260px, 1fr) minmax(260px, 1fr);
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
      <div class="stat"><b>{len(entries):,}</b><span>rule entries</span></div>
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
