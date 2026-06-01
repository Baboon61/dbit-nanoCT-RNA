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


def format_metric(key, value):
    if key.endswith("_bytes") or key == "bytes":
        return format_bytes(value)
    return format_value(value)


def normalize_rule(rule):
    if rule.startswith("matrix_bin_"):
        return "create_matrix_bins"
    return RULE_ALIASES.get(rule, rule)


def entry_context_label(entry):
    context = entry.get("context") or {}
    ordered_keys = ["sample", "modality", "barcode", "number", "lane", "suffix", "ext", "matrix"]
    parts = []
    for key in ordered_keys:
        if key in context:
            parts.append(f"{key}: {context[key]}")
    for key in sorted(set(context) - set(ordered_keys)):
        parts.append(f"{key}: {context[key]}")
    return "; ".join(parts) or "global"


def render_key_values(values, class_name="kv"):
    if not values:
        return '<span class="muted">none</span>'
    rows = []
    for key in sorted(values):
        rows.append(
            "<tr>"
            f"<th>{human_label(key)}</th>"
            f"<td>{format_metric(key, values[key])}</td>"
            "</tr>"
        )
    return f'<table class="{class_name}"><tbody>' + "".join(rows) + "</tbody></table>"


def render_outputs(outputs):
    if not outputs:
        return '<span class="muted">none</span>'
    items = "".join(f"<li>{escape(path)}</li>" for path in outputs)
    return f"<ul>{items}</ul>"


def render_entry(entry):
    target = entry.get("target_label") or entry.get("rule")
    context = entry_context_label(entry)
    metrics = render_key_values(entry.get("metrics") or {})
    runtime = render_key_values(entry.get("runtime") or {})
    outputs = render_outputs(entry.get("outputs") or [])
    return f"""
      <article class="entry">
        <div class="entry-head">
          <h3>{escape(target)}</h3>
          <p>{escape(context)}</p>
        </div>
        <div class="entry-grid">
          <section>
            <h4>Metrics</h4>
            {metrics}
          </section>
          <section>
            <h4>Runtime</h4>
            {runtime}
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


def build_html(report, rule_order):
    entries = report.get("rules") or []
    order_index = {rule: index for index, rule in enumerate(rule_order)}
    groups = defaultdict(list)
    for entry in entries:
        groups[normalize_rule(entry.get("rule", "unknown"))].append(entry)

    total_outputs = sum(len(entry.get("outputs") or []) for entry in entries)
    metric_entries = sum(1 for entry in entries if entry.get("metrics"))
    workflow = report.get("workflow", "workflow")

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
              {''.join(render_entry(entry) for entry in rule_entries)}
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
    .entry-head h3 {{
      font-size: 1rem;
      font-weight: 650;
    }}
    .entry-head p {{
      color: var(--muted);
      font-size: 0.9rem;
      text-align: right;
      overflow-wrap: anywhere;
    }}
    .entry-grid {{
      display: grid;
      grid-template-columns: minmax(220px, 1.2fr) minmax(180px, 0.8fr) minmax(260px, 1fr);
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
    .muted {{
      color: var(--muted);
      font-size: 0.9rem;
    }}
    @media (max-width: 860px) {{
      .entry-head {{
        display: block;
      }}
      .entry-head p {{
        text-align: left;
        margin-top: 4px;
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
    <p class="subtitle">Rule metrics from {escape(report.get("processedData_dir", ""))}, grouped and ordered by the workflow rule sequence.</p>
  </header>
  <main>
    <section class="summary">
      <div class="stat"><b>{len(entries):,}</b><span>rule entries</span></div>
      <div class="stat"><b>{len(groups):,}</b><span>rules represented</span></div>
      <div class="stat"><b>{metric_entries:,}</b><span>entries with metrics</span></div>
      <div class="stat"><b>{total_outputs:,}</b><span>tracked outputs</span></div>
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
