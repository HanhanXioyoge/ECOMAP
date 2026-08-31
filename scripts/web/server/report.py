"""
HTML report renderer for /api/jobs/{jid}/report (B14 / T30).

Kept deliberately tiny — Jinja2 template with a single ``render_report``
entry point. The route handler in ``routes_jobs_report`` is responsible
for fetching the job record and feeding it to ``render_report``.
"""
from __future__ import annotations

import json
from typing import Any

from jinja2 import Template


_REPORT_TPL = Template(
    """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ECOMAP Job {{ jid }}</title>
<style>
body { font-family: Inter, system-ui, sans-serif; max-width: 900px; margin: 24px auto; padding: 0 16px; color: #1A1A1A; }
h1 { font-size: 20px; border-bottom: 1px solid #E5E7EB; padding-bottom: 8px; }
h2 { font-size: 16px; margin-top: 24px; }
pre { background: #F9FAFB; padding: 12px; border: 1px solid #E5E7EB; overflow-x: auto; font-size: 12px; }
table { border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 14px; }
th, td { border: 1px solid #E5E7EB; padding: 6px 10px; text-align: left; }
th { background: #F9FAFB; font-weight: 500; }
.meta dt { font-weight: 600; display: inline-block; width: 110px; color: #4B5563; }
.meta dd { display: inline; margin: 0; }
.meta div { margin: 2px 0; }
</style>
</head>
<body>
<h1>ECOMAP Job Report</h1>
<dl class="meta">
  <div><dt>Job ID:</dt><dd>{{ jid }}</dd></div>
  <div><dt>Status:</dt><dd>{{ status }}</dd></div>
  <div><dt>Track:</dt><dd>{{ track }}</dd></div>
  <div><dt>Started:</dt><dd>{{ started_at or "-" }}</dd></div>
  <div><dt>Ended:</dt><dd>{{ ended_at or "-" }}</dd></div>
</dl>

<h2>Log</h2>
<pre>{{ log_text }}</pre>

<h2>Result</h2>
{% if result_table %}
<table>
  <thead><tr>{% for c in result_columns %}<th>{{ c }}</th>{% endfor %}</tr></thead>
  <tbody>
  {% for row in result_table %}
    <tr>{% for v in row %}<td>{{ v }}</td>{% endfor %}</tr>
  {% endfor %}
  </tbody>
</table>
{% else %}
<pre>{{ result_json }}</pre>
{% endif %}
</body>
</html>
"""
)


def render_report(
    jid: str,
    status: str,
    track: str,
    started_at: str,
    ended_at: str,
    log: Any,
    result: Any,
) -> str:
    """Render the job's report to a self-contained HTML string.

    ``log`` may be a string or an iterable of log lines; ``result`` is
    either a tabular ``{"columns": [...], "rows": [...]}`` payload or any
    other dict/list that gets JSON-serialised into a ``<pre>`` block.
    """
    if isinstance(log, str):
        log_text = log
    else:
        log_text = "\n".join(str(line) for line in (log or []))

    result_table: list | None = None
    result_columns: list[str] = []
    result_json = ""
    if isinstance(result, dict) and "columns" in result and "rows" in result:
        result_columns = list(result.get("columns") or [])
        result_table = list(result.get("rows") or [])
    else:
        result_json = json.dumps(result, indent=2, default=str)

    return _REPORT_TPL.render(
        jid=jid,
        status=status,
        track=track,
        started_at=started_at,
        ended_at=ended_at,
        log_text=log_text,
        result_table=result_table,
        result_columns=result_columns,
        result_json=result_json,
    )
