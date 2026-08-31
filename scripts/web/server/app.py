# server/app.py
"""FastAPI application for the MDP Web UI backend.

Routes (v1):
    GET  /api/i18n/{lang}            -> i18n string table for {zh, en}
    POST /api/model                  -> upload a model file, return model_info + modelId
    GET  /api/algorithms             -> algorithm registry filtered by ?modelType=...
    POST /api/jobs                   -> start a background job (FSEOF, etc.)
    GET  /api/jobs/{jid}             -> job status snapshot
    GET  /api/jobs/{jid}/events      -> SSE event stream for one job
    GET  /api/jobs/{jid}/candidates  -> final FSEOF {columns, rows, log}
    GET  /api/jobs/{jid}/export.csv  -> CSV download of the candidates table

The static frontend (web/, filled by Task 5) is mounted LAST so it never
shadows the /api/* routes.
"""
import asyncio
import csv
import io
import json
import os
import shutil
import uuid
from pathlib import Path
from queue import Empty as _QueueEmpty

from fastapi import FastAPI, UploadFile, File, HTTPException, Request
from fastapi.responses import StreamingResponse, Response, FileResponse
from fastapi.staticfiles import StaticFiles

import jobs
import matlab_bridge as mb
import projects as project_store
from routes_algorithms import algorithms_router
from routes_jobs_artifacts import jobs_artifacts_router
from routes_jobs_report import jobs_report_router
from routes_recon_init import recon_init_router
from routes_recon_param import recon_param_router

# Repo root = parent of scripts/ (this file lives at scripts/web/server/app.py;
# four dirname() hops climb server -> web -> scripts -> <repo root>).
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
UPLOADS = os.path.join(os.path.dirname(__file__), "_uploads")
os.makedirs(UPLOADS, exist_ok=True)

# In-process registry: modelId -> stored file path. Lives for the process
# lifetime only (v1 is single-user localhost).
_models: dict[str, str] = {}
_uploads: dict[str, str] = {}

app = FastAPI(title="ECOMAP")


@app.get("/api/i18n/{lang}")
def get_i18n(lang: str):
    if lang not in ("zh", "en"):
        raise HTTPException(400, "bad lang")
    path = os.path.join(ROOT, "scripts", "web", "matlab", "i18n", f"{lang}.json")
    with open(path, encoding="utf-8") as f:  # UTF-8 per global constraints.
        return json.load(f)


@app.post("/api/projects")
def create_project(payload: dict):
    try:
        return project_store.create_project(payload)
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc


@app.get("/api/projects")
def list_projects(limit: int = 20):
    n = max(0, min(int(limit), 100))
    return {"projects": project_store.list_projects(n)}


@app.get("/api/projects/recent")
def recent_projects(limit: int = 3):
    n = max(0, min(int(limit), 20))
    return {"projects": project_store.list_projects(n)}


@app.get("/api/projects/{project_id}")
def get_project(project_id: str):
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(404, "unknown project")
    return project


@app.get("/api/projects/{project_id}/files")
def get_project_files(project_id: str):
    try:
        return project_store.get_project_files(project_id)
    except KeyError as exc:
        raise HTTPException(404, "unknown project") from exc


@app.get("/api/projects/{project_id}/reconstruction/state")
def get_project_reconstruction_state(project_id: str):
    try:
        return project_store.get_reconstruction_state(project_id)
    except KeyError as exc:
        raise HTTPException(404, "unknown project") from exc


@app.get("/api/projects/{project_id}/file/{relative_path:path}")
def get_project_file(project_id: str, relative_path: str):
    try:
        path = project_store.resolve_project_file(project_id, relative_path)
    except KeyError as exc:
        raise HTTPException(404, "unknown project file") from exc
    return FileResponse(path)


@app.post("/api/projects/{project_id}/models")
async def upload_project_model(project_id: str, file: UploadFile = File(...)):
    try:
        path = project_store.project_model_upload_path(project_id, file.filename or "model.mat")
    except KeyError as exc:
        raise HTTPException(404, "unknown project") from exc
    with open(path, "wb") as out:
        shutil.copyfileobj(file.file, out)
    mid = uuid.uuid4().hex[:12]
    _models[mid] = str(path)
    rel = Path("models") / path.name
    return {
        "modelId": mid,
        "uploadId": mid,
        "name": path.name,
        "path": str(path),
        "projectRelativePath": rel.as_posix(),
    }


@app.put("/api/projects/{project_id}/params")
def update_project_params(project_id: str, payload: dict):
    try:
        project = project_store.update_project_params(project_id, payload.get("params") or payload)
        if payload.get("loadManager", True):
            manager_path = project_store.parameter_manager_path(project_id)
            project["parameterManagerLoaded"] = mb.load_parameter_manager(manager_path)
        return project
    except KeyError as exc:
        raise HTTPException(404, "unknown project") from exc
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, str(exc)) from exc


@app.post("/api/model")
async def post_model(file: UploadFile = File(...)):
    # Unique subdirectory per upload keeps original filenames collision-free.
    # Strip path components to prevent traversal (e.g. "../../etc/...").
    mid = uuid.uuid4().hex[:12]
    d = os.path.join(UPLOADS, mid)
    os.makedirs(d, exist_ok=True)
    safe_name = os.path.basename(file.filename or "model.mat")
    path = os.path.join(d, safe_name)
    with open(path, "wb") as out:
        shutil.copyfileobj(file.file, out)
    _models[mid] = path
    info = mb.model_info(path)
    return {"modelId": mid, **info}


# Generic upload endpoint used by recon.js (model files) and calib.js
# (proteomics data). Returns {uploadId, path, name, modelId?} so both call
# sites can keep their existing field names: recon reads up.modelId, calib
# reads up.uploadId. modelId is only populated when the upload looks like a
# metabolic model (so the recon wizard can advance straight to convert).
_MODEL_EXT = {".mat", ".xml", ".json", ".yml", ".yaml"}


@app.post("/api/uploads")
async def post_upload(file: UploadFile = File(...)):
    """Store an uploaded file and return its handle.

    Works for both model files and arbitrary supporting data. ``modelId`` is
    only populated for recognised model extensions and is best-effort: the
    underlying ``mdpModelInfo`` bridge can still fail (e.g. when a Parameter
    Manager context is required). On bridge failure the upload itself is
    preserved so the caller can retry or attach it to a project first.
    """
    uid = uuid.uuid4().hex[:12]
    d = os.path.join(UPLOADS, uid)
    os.makedirs(d, exist_ok=True)
    safe_name = os.path.basename(file.filename or "upload.bin")
    path = os.path.join(d, safe_name)
    with open(path, "wb") as out:
        shutil.copyfileobj(file.file, out)
    resp = {"uploadId": uid, "path": path, "name": safe_name}
    _uploads[uid] = path
    if Path(safe_name).suffix.lower() in _MODEL_EXT:
        _models[uid] = path
        resp["modelId"] = uid
        try:
            info = mb.model_info(path)
            resp.update(info)
        except Exception as e:  # noqa: BLE001 - bridge failure must not break upload
            resp["model_info_error"] = str(e)
    return resp


def _require(payload: dict, key: str):
    value = payload.get(key)
    if value in (None, "", [], {}):
        raise HTTPException(400, f"missing field: {key}")
    return value


def _data_value(data: dict, method: str, field: str):
    value = ((data or {}).get(method) or {}).get(field)
    if value in (None, ""):
        raise HTTPException(400, f"missing data: {method}.{field}")
    return _uploads.get(value, value)


def _project_model_path(payload: dict, params: dict, model_id: str) -> str | None:
    project_id = payload.get("projectId")
    if not project_id:
        return None
    candidates = []
    if model_id:
        candidates.append(model_id)
        if not str(model_id).startswith("models/"):
            candidates.append(f"models/{os.path.basename(str(model_id))}")
    initial = params.get("InitialModel")
    if initial:
        candidates.append(str(initial))
        if not str(initial).startswith("models/"):
            candidates.append(f"models/{os.path.basename(str(initial))}")
    for candidate in candidates:
        try:
            return str(project_store.resolve_project_file(str(project_id), str(candidate)))
        except KeyError:
            continue
    return None


def _ensure_project_parameter_manager(payload: dict) -> str | None:
    project_id = payload.get("projectId")
    if not project_id:
        return None
    manager_path = project_store.parameter_manager_path(str(project_id))
    mb.load_parameter_manager(manager_path)
    return str(manager_path)


def _with_web_model_id(result: dict) -> dict:
    """Expose MATLAB's ``model_id`` under the frontend's ``modelId`` name."""
    normalized = dict(result)
    if normalized.get("model_id") and not normalized.get("modelId"):
        normalized["modelId"] = normalized["model_id"]
    return normalized


def _with_web_ec_model_ids(result: dict, topology: str | None = None) -> dict:
    """Expose MATLAB's ecModel handles under the frontend's camelCase names."""
    normalized = dict(result)
    ec_ids = normalized.get("ecModelIds") or normalized.get("ec_model_ids") or {}
    if ec_ids and not normalized.get("ecModelIds"):
        normalized["ecModelIds"] = ec_ids
    if ec_ids and not normalized.get("ecModelId"):
        if topology and topology in ec_ids:
            normalized["ecModelId"] = ec_ids[topology]
        else:
            normalized["ecModelId"] = next(iter(ec_ids.values()), "")
    return normalized


def _with_web_merged_ec_model_ids(result: dict) -> dict:
    normalized = dict(result)
    ec_ids = normalized.get("ecModelIds") or normalized.get("merged_ec_model_ids") or {}
    if ec_ids and not normalized.get("ecModelIds"):
        normalized["ecModelIds"] = ec_ids
    if ec_ids and not normalized.get("ecModelId"):
        normalized["ecModelId"] = next(iter(ec_ids.values()), "")
    return normalized


def _selected_topologies(params: dict) -> list[str]:
    allowed = {"basic", "isozyme", "integrated"}
    raw = params.get("topologies")
    if isinstance(raw, list):
        candidates = [str(item) for item in raw]
    else:
        topology = str(params.get("topology", "integrated"))
        candidates = ["basic", "isozyme", "integrated"] if topology == "all" else [topology]
    selected = []
    for topology in candidates:
        if topology not in allowed:
            raise HTTPException(400, f"invalid topology: {topology}")
        if topology not in selected:
            selected.append(topology)
    if not selected:
        raise HTTPException(400, "at least one topology is required")
    return selected


def _annotation_options(params: dict) -> dict:
    raw = params.get("annotationOptions") if isinstance(params.get("annotationOptions"), dict) else {}
    stages = params.get("annotationStages") or ["ec", "metabolite"]
    sources = raw.get("metaboliteSources", params.get("metaboliteSources", ["A", "B", "C"]))
    if not isinstance(sources, list):
        sources = [sources]
    selected_sources = []
    for source in [str(item).upper() for item in sources if str(item)]:
        if source not in {"A", "B", "C"}:
            raise HTTPException(400, f"invalid metabolite source: {source}")
        if source not in selected_sources:
            selected_sources.append(source)
    return {
        "runComplexAnnotation": bool(raw.get(
            "runComplexAnnotation",
            params.get("runComplexAnnotation", "complex" in stages),
        )),
        "runEcAnnotation": bool(raw.get(
            "runEcAnnotation",
            params.get("runEcAnnotation", "ec" in stages),
        )),
        "metaboliteSources": selected_sources,
        "runMetaNetXIntegration": bool(raw.get(
            "runMetaNetXIntegration",
            params.get("runMetaNetXIntegration", True),
        )),
    }


def _ordered_ec_model_ids(payload: dict) -> list[str]:
    ec_model_ids = payload.get("ecModelIds") or {}
    if isinstance(ec_model_ids, dict) and ec_model_ids:
        ordered = []
        for topology in ("integrated", "basic", "isozyme"):
            value = ec_model_ids.get(topology)
            if value:
                ordered.append(value)
        ordered.extend(value for key, value in ec_model_ids.items() if key not in {"integrated", "basic", "isozyme"} and value)
        return ordered
    return [_require(payload, "ecModelId")]


def _ec_model_ids_by_topology(payload: dict) -> dict[str, str]:
    raw = payload.get("ecModelIds") or {}
    if not isinstance(raw, dict):
        return {}
    return {
        str(topology): str(model_id)
        for topology, model_id in raw.items()
        if topology in {"integrated", "isozyme", "basic"} and model_id
    }


def _kcat_reference_topology(params: dict, ec_model_ids: dict[str, str]) -> str:
    requested = params.get("kcatReferenceTopology") or params.get("customKcatRxnNameType")
    if requested and requested in ec_model_ids:
        return str(requested)
    for topology in ("integrated", "isozyme", "basic"):
        if topology in ec_model_ids:
            return topology
    return str(requested or "")


def _kcat_ec_model_id(payload: dict, params: dict) -> str:
    ec_model_ids = _ec_model_ids_by_topology(payload)
    reference = _kcat_reference_topology(params, ec_model_ids)
    if reference and reference in ec_model_ids:
        return ec_model_ids[reference]
    return _require(payload, "ecModelId")


def _kcat_prediction_model(params: dict) -> str:
    models = params.get("deepLearningModels") or ["DLKcat", "UniKP", "CatPred"]
    requested = params.get("kcatPredictionModel")
    if requested and requested in models:
        return str(requested)
    if "CatPred" in models:
        return "CatPred"
    return str(models[-1])


def _kcat_merge_options(params: dict, ec_model_ids: dict[str, str]) -> dict:
    reference = _kcat_reference_topology(params, ec_model_ids)
    return {
        "kcatReferenceTopology": reference,
        "kcatPredictionModel": _kcat_prediction_model(params),
        "customKcatRxnNameType": params.get("customKcatRxnNameType") or reference,
        "medianThreshold": params.get("medianThreshold", 5),
        "useLoggedMedian": bool(params.get("useLoggedMedian", True)),
    }


@app.post("/api/recon/run")
def recon_run(payload: dict):
    """Run one Reconstruction action through the MATLAB bridge."""
    action = _require(payload, "action")
    params = payload.get("params") or {}
    try:
        if action == "init":
            project_name = params.get("projectName") or payload.get("projectId") or "ECOMAP_Project"
            project_path = params.get("projectPath") or os.path.join(UPLOADS, "projects", str(project_name))
            return {"ok": True, "result": mb.init_project(str(project_name), str(project_path))}
        manager_path = _ensure_project_parameter_manager(payload)
        if action == "load":
            model_id = _require(payload, "modelId")
            path = _models.get(model_id)
            if not path:
                path = _project_model_path(payload, params, str(model_id))
            if not path:
                raise HTTPException(404, "unknown modelId")
            result = mb.load_model(path, params.get("modeltype", "Tradition"), manager_path)
            return {"ok": True, "result": _with_web_model_id(result)}
        if action == "convert":
            model_id = _require(payload, "modelId")
            topologies = _selected_topologies(params)
            if len(topologies) > 1:
                results = {}
                ec_model_ids = {}
                for topo in topologies:
                    result = _with_web_ec_model_ids(mb.convertec_model(model_id, topo, manager_path), topo)
                    results[topo] = result
                    if result.get("ecModelId"):
                        ec_model_ids[topo] = result["ecModelId"]
                return {"ok": True, "result": {"topologies": results, "ecModelIds": ec_model_ids}}
            topology = topologies[0]
            result = mb.convertec_model(model_id, topology, manager_path)
            return {"ok": True, "result": _with_web_ec_model_ids(result, topology)}
        if action == "annotate":
            ec_ids = _ordered_ec_model_ids(payload)
            stages = params.get("annotationStages") or ["ec", "metabolite"]
            return {"ok": True, "result": mb.annotate(ec_ids, stages, manager_path, _annotation_options(params))}
        if action == "deepLearningKcat":
            ec_model_id = _kcat_ec_model_id(payload, params)
            models = params.get("deepLearningModels") or ["DLKcat", "UniKP", "CatPred"]
            return {"ok": True, "result": mb.dl_predict(ec_model_id, models, manager_path)}
        if action == "compare":
            ec_model_id = _kcat_ec_model_id(payload, params)
            models = params.get("deepLearningModels") or ["DLKcat", "UniKP", "CatPred"]
            return {"ok": True, "result": mb.kcat_compare(ec_model_id, models, params.get("complexNames") or [], manager_path)}
        if action == "merge":
            ec_model_ids = _ec_model_ids_by_topology(payload)
            ec_ids = _ordered_ec_model_ids(payload)
            result = mb.kcat_merge(ec_ids, bool(params.get("useCustomKcatFile", False)), _kcat_merge_options(params, ec_model_ids), manager_path)
            return {"ok": True, "result": _with_web_merged_ec_model_ids(result)}
        if action == "growth":
            ec_model_id = _require(payload, "ecModelId")
            return {"ok": True, "result": mb.growth_predict(ec_model_id, params.get("carbonSource", "glucose"), params.get("biomassReaction", "biomass"), manager_path)}
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, str(exc)) from exc
    raise HTTPException(400, f"unknown reconstruction action: {action}")


@app.post("/api/calib/run")
def calib_run(payload: dict):
    """Run selected Calibration methods in the caller-supplied fixed order."""
    ec_model_id = _require(payload, "ecModelId")
    methods = payload.get("methods") or []
    params = payload.get("params") or {}
    data = payload.get("data") or {}
    if not methods:
        raise HTTPException(400, "missing field: methods")

    results = []
    try:
        for method in methods:
            if method == "sluice":
                ex_rxns = params.get("exRxns") or []
                if isinstance(ex_rxns, str):
                    ex_rxns = [x.strip() for x in ex_rxns.split(",") if x.strip()]
                results.append({"method": method, "result": mb.apply_sluice(ec_model_id, ex_rxns)})
            elif method == "kcatRepo":
                results.append({"method": method, "result": mb.kcat_repo_init(ec_model_id)})
            elif method == "sensitivity":
                _data_value(data, "sensitivity", "glucoseUptake")
                _data_value(data, "sensitivity", "growthRate")
                results.append({
                    "method": method,
                    "result": mb.sensitivity_tuning(
                        ec_model_id,
                        params.get("glucoseExchange", "EX_glc__D_e"),
                        float(params.get("targetGrowth", 0.0)),
                        float(params.get("factor", 1.0)),
                        bool(params.get("multiCondition", False)),
                    ),
                })
            elif method == "bayesian":
                scenario = _data_value(data, "bayesian", "scenarioData")
                results.append({
                    "method": method,
                    "result": mb.bayesian(
                        ec_model_id,
                        [scenario],
                        int(params.get("maxIter", 200)),
                        int(params.get("proc", 4)),
                        int(params.get("numPerGen", 100)),
                        float(params.get("rejectNum", 0.05)),
                        bool(params.get("runGauksAfterBayesian", False)),
                    ),
                })
            elif method == "presto":
                files = {
                    "proteomics": _data_value(data, "presto", "proteomics"),
                    "growth": _data_value(data, "presto", "growth"),
                    "total_protein": _data_value(data, "presto", "totalProtein"),
                }
                results.append({"method": method, "result": mb.presto(ec_model_id, files)})
            elif method == "gauks":
                _data_value(data, "gauks", "unconstrainedMaxGrowth")
                gem_model_id = _require(payload, "gemModelId")
                results.append({
                    "method": method,
                    "result": mb.gauks(ec_model_id, gem_model_id, params.get("biomassReaction", "biomass")),
                })
            else:
                raise HTTPException(400, f"unknown calibration method: {method}")
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, str(exc)) from exc

    return {"ok": True, "results": results}


# /api/algorithms (T31) lives in routes_algorithms.py — see algorithms_router.


# ---------------------------------------------------------------------------
# Jobs (async FSEOF runs + SSE + candidates + CSV export)
# ---------------------------------------------------------------------------


@app.post("/api/jobs")
def post_job(payload: dict):
    mid = payload.get("modelId")
    if not mid or mid not in _models:
        raise HTTPException(404, "unknown modelId")
    algos = payload.get("algos") or []
    spec = None
    if payload.get("algo"):
        spec = {"id": payload.get("algo"), "params": payload.get("params") or {}}
    elif algos:
        spec = algos[0]
    if not spec or not spec.get("id"):
        raise HTTPException(400, "msg_select_algo")
    algo = str(spec["id"]).lower()
    if algo not in {"fseof", "optknock", "optforce", "oko", "okoplus", "okoplus-build"}:
        raise HTTPException(400, "unsupported algorithm")
    p = dict(spec.get("params") or {})
    if algo == "okoplus" and not p.get("interval_path"):
        source_jid = p.get("intervalJobId") or p.get("interval_job_id")
        predictor = str(p.get("predictor", "UniKP")).lower()
        if source_jid:
            interval_path = jobs.get_workspace_file(
                str(source_jid), f"intervals_{predictor}.csv")
            if not interval_path:
                raise HTTPException(400, "interval artifact not found")
            p["interval_path"] = interval_path
    project_id = payload.get("projectId")
    if project_id and algo == "okoplus-build":
        p["manager_path"] = str(project_store.parameter_manager_path(str(project_id)))

    model_ref = _models[mid]
    # Upload handles store a path, whereas design bridges consume a MATLAB
    # registry id.  Register on first design submission and cache the handle.
    if os.path.isfile(str(model_ref)):
        try:
            loaded = mb.load_model(str(model_ref), str(payload.get("modelType", "Tradition")))
            model_ref = loaded.get("model_id") or loaded.get("modelId")
            if not model_ref:
                raise RuntimeError("MATLAB did not return model_id")
            _models[mid] = model_ref
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(500, str(exc)) from exc
    # MATLAB errors (err_no_target, err_no_biomass, Gurobi errors, ...) surface
    # asynchronously via SSE /api/jobs/{jid}/events and the status endpoint;
    # jobs.create_job only registers state and spawns a daemon thread, so the
    # synchronous POST /api/jobs always succeeds (HTTP 200) as long as the input
    # shape above is OK.
    jid = jobs.create_job(
        model_ref,
        str(payload.get("biomass", "")),
        str(payload.get("target", "")),
        algo,
        p,
    )
    return {"jobId": jid}


@app.get("/api/jobs/{jid}")
def job_status(jid: str):
    j = jobs.get_job(jid)
    if not j:
        raise HTTPException(404, "unknown job")
    return {"status": j["status"], "error": j["error"]}


@app.get("/api/jobs/{jid}/result")
def job_result(jid: str):
    j = jobs.get_job(jid)
    if not j:
        raise HTTPException(404, "unknown job")
    if j["status"] != "done":
        raise HTTPException(409, "job not complete")
    return j.get("result") or {}


@app.get("/api/jobs")
def list_jobs(limit: int = 3):
    """Return the most recent N jobs for the home-view recent works card.

    Limit is clamped to [0, 20] to keep the response small.
    """
    n = max(0, min(int(limit), 20))
    return {"jobs": jobs.list_jobs(n)}


@app.get("/api/jobs/{jid}/events")
async def job_events(jid: str, request: Request):
    q = jobs.get_queue(jid)
    if q is None:
        raise HTTPException(404, "unknown job")

    async def gen():
        while True:
            if await request.is_disconnected():
                break
            try:
                ev = await asyncio.to_thread(q.get, timeout=0.5)
            except _QueueEmpty:  # poll timeout -> try again
                continue
            yield f"data: {json.dumps(ev)}\n\n"
            if ev.get("type") in ("done", "error"):
                break

    return StreamingResponse(gen(), media_type="text/event-stream")


def _reshape_rows(res: dict) -> list:
    """Normalise the MATLAB rows payload into a list-of-lists, one per row.

    The MATLAB adapter serialises a 13x7 `table` via `table2cell` then
    `jsonencode`. Depending on MATLAB's encoder version the result may be:
      * a nested list-of-lists already in row-major form (13 rows, each 7 cols)
      * a flat list of length 13*7 = 91 that needs reshaping
      * a column-major flat list where 7 entries per "column" appear back-to-back,
        producing a layout that reads as 7 rows of 13 cols when naively reshaped

    Accept any of these and emit the canonical 13-row × 7-col shape.
    """
    cols = res.get("columns") or []
    rows = res.get("rows") or []
    n = len(cols)
    if n == 0 or not rows:
        return []

    # Already nested list-of-lists? Inspect its actual width.
    if isinstance(rows[0], list):
        width = len(rows[0])
        if width == n:
            return [(r + [None] * n)[:n] for r in rows]
        # Column-major nested: width is the number of records (13) but each
        # "row" actually represents one column across all fields.
        if len(rows) == n:
            return [[row[i] for row in rows] for i in range(width)]
        # Fallback: trust whatever was given, trim/pad defensively.
        return [(r + [None] * n)[:n] for r in rows]

    # Flat list — could be row-major (n*num_records) or column-major (num_records*n).
    total = len(rows)
    if total == 0:
        return []
    if total % n == 0:
        # Row-major: chunks of n.
        num = total // n
        return [rows[i:i + n] for i in range(0, total, n)]
    # Try the column-major shape: total = num_records * n_records.
    # We don't know num_records; assume width = n (fields) and length = num_records.
    # If num_records is reached, num_records * n == total — so width must be n
    # but rows had been encoded in column order, giving width = n and the
    # wrong row count. We can't recover without metadata; fall through to a
    # best-effort split into num_records=n groups of size n.
    if total % n == 0:
        num = total // n
        return [rows[i:i + n] for i in range(0, total, n)]
    # Last resort: truncate.
    m = (total // n) * n
    return [rows[i:i + n] for i in range(0, m, n)]


@app.get("/api/jobs/{jid}/candidates")
def job_candidates(jid: str):
    res = jobs.get_result(jid)
    if res is None:
        raise HTTPException(404, "no result")
    shaped = dict(res)
    shaped["rows"] = _reshape_rows(res)
    return shaped


@app.get("/api/jobs/{jid}/export.csv")
def job_export(jid: str):
    res = jobs.get_result(jid)
    if res is None:
        raise HTTPException(404, "no result")
    shaped_rows = _reshape_rows(res)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(res["columns"])
    w.writerows(shaped_rows)
    return Response(
        buf.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=candidates.csv"},
    )


# ---------------------------------------------------------------------------
# Routers added by B14 / T27-T31 (recon, jobs artifacts + report, algorithms)
# ---------------------------------------------------------------------------
app.include_router(recon_init_router)
app.include_router(recon_param_router)
app.include_router(jobs_artifacts_router)
app.include_router(jobs_report_router)
app.include_router(algorithms_router)


# Static frontend mounted LAST so it never shadows /api/*. Conditional:
# only mount when the directory exists so the app boots before the frontend lands.
_WEB_DIR = os.path.join(ROOT, "scripts", "web")
if os.path.isdir(_WEB_DIR) and os.path.isfile(os.path.join(_WEB_DIR, "index.html")):
    app.mount("/", StaticFiles(directory=_WEB_DIR, html=True), name="web")


# NOTE: host/port are fixed for v1 (uvicorn server.app:app --host 127.0.0.1
# --port 8000). MDP_HOST / configurable port to be wired in a later task.
