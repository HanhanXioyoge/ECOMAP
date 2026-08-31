"""
Jobs router — T30: GET /api/jobs/{jid}/report.

Renders the per-job HTML report via ``report.render_report``. The job
record is fetched from ``jobs.get_job``; missing keys (e.g. ``track``,
``started_at``, ``log``) fall back to empty strings so the report
template always renders cleanly.
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException
from fastapi.responses import Response

import jobs
from report import render_report

jobs_report_router = APIRouter()


@jobs_report_router.get("/api/jobs/{jid}/report")
def job_report(jid: str):
    """Render an HTML report for one job."""
    j = jobs.get_job(jid)
    if j is None:
        raise HTTPException(404, "unknown job")
    html = render_report(
        jid=jid,
        status=j.get("status", "?"),
        track=j.get("track", "?"),
        started_at=j.get("started_at", ""),
        ended_at=j.get("ended_at", ""),
        log=j.get("log", ""),
        result=j.get("result", {}),
    )
    return Response(html, media_type="text/html")
