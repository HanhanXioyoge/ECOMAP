"""
Recon track router — T27: POST /api/recon/init.

Initialises a new ECOMAP project and returns its id + default parameter
template. The MATLAB ``mdpInitProject`` call is delegated to
``matlab_bridge.init_project``, which already unwraps the bridge envelope
and raises ``BridgeContractError`` on failure.
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException

import matlab_bridge as mb

recon_init_router = APIRouter()


@recon_init_router.post("/api/recon/init")
def recon_init(payload: dict):
    """Initialise a new ECOMAP project; return its id + default params."""
    try:
        project_name = payload["project_name"]
        project_path = payload["project_path"]
    except KeyError as exc:
        raise HTTPException(400, f"missing field: {exc.args[0]}") from exc
    return mb.init_project(project_name, project_path)
