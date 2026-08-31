"""
Recon track router — T28: GET + PUT /api/recon/param/{projectId}.

Reads and writes the per-project parameter JSON used to drive the
ParameterManagement.m MATLAB script. The store is a tiny on-disk
key/value map: one file per projectId under ``scripts/web/server/.uploads/``.
"""
from __future__ import annotations

import json
import os

from fastapi import APIRouter, HTTPException

# Per-project parameter JSON files live under this directory. Leading dot
# keeps the directory out of any "uploads" UI listings. Created at import
# time so PUT can write into it without an extra mkdir round-trip.
PARAMS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".uploads")
os.makedirs(PARAMS_DIR, exist_ok=True)


def _param_path(project_id: str) -> str:
    """Absolute path to the JSON file holding the saved params for a project."""
    # ``os.path.basename`` strips any path components the caller might
    # inject via the URL — we never want to escape PARAMS_DIR.
    safe = os.path.basename(project_id)
    return os.path.join(PARAMS_DIR, f"{safe}.params.json")


recon_param_router = APIRouter()


@recon_param_router.get("/api/recon/param/{project_id}")
def recon_param_get(project_id: str):
    """Return the previously-saved parameter JSON for one project."""
    path = _param_path(project_id)
    if not os.path.exists(path):
        raise HTTPException(404, f"no params for project_id={project_id!r}")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


@recon_param_router.put("/api/recon/param/{project_id}")
def recon_param_put(project_id: str, payload: dict):
    """Persist the parameter JSON for one project to disk."""
    path = _param_path(project_id)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return {"ok": True, "project_id": project_id}
