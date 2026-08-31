"""
Jobs router — T29: POST /api/jobs/{jid}/cancel + GET .../ecModel + GET .../kcatRepo.

Provides:

- ``cancel`` — flips the in-memory job status to ``cancelled``; returns 200
  on success and 404 on an unknown jid.
- ``ecModel`` and ``kcatRepo`` — stream the named artifact back to the
  client via ``FileResponse``. Both rely on the per-job ``artifacts`` map
  populated by the worker thread.
"""
from __future__ import annotations

import os

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

import jobs

jobs_artifacts_router = APIRouter()


@jobs_artifacts_router.post("/api/jobs/{jid}/cancel")
def job_cancel(jid: str):
    """Flip a running/queued job to ``cancelled``."""
    if not jobs.cancel_job(jid):
        raise HTTPException(404, "unknown job")
    return {"jid": jid, "status": "cancelled"}


def _serve_artifact(jid: str, name: str):
    """Return a ``FileResponse`` for the named artifact, or raise 404."""
    path = jobs.get_workspace_file(jid, name)
    if path is None or not os.path.exists(path):
        raise HTTPException(404, f"no {name} for jid={jid!r}")
    return FileResponse(path, filename=os.path.basename(path))


@jobs_artifacts_router.get("/api/jobs/{jid}/ecModel")
def job_ec_model(jid: str):
    """Serve the ecModel .mat artifact for a job."""
    return _serve_artifact(jid, "ecModel")


@jobs_artifacts_router.get("/api/jobs/{jid}/kcatRepo")
def job_kcat_repo(jid: str):
    """Serve the KcatRepo .mat artifact for a job."""
    return _serve_artifact(jid, "kcatRepo")


@jobs_artifacts_router.get("/api/jobs/{jid}/artifacts/{name}")
def job_named_artifact(jid: str, name: str):
    """Serve one explicitly registered job artifact by its stable name."""
    if "/" in name or "\\" in name or name in {".", "..", "run_dir"}:
        raise HTTPException(400, "invalid artifact name")
    return _serve_artifact(jid, name)
