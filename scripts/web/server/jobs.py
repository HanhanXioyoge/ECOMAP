"""In-process job subsystem for the ECOMAP Web UI.

Jobs are tracked in two module-level dicts:

    _jobs   jobId -> dict(status, result, error, algo, params, ...)
    _queues jobId -> queue.Queue of event dicts

A short-lived background thread runs the actual MATLAB call under
``_worker_lock`` so the single MATLAB engine is never entered concurrently.
Events (log lines, terminal ``done``/``error``) are pushed into the per-job
queue for SSE consumers to drain.

Per the v1 constraints: in-memory only, single-user, localhost. No
persistence, no auth, no cleanup daemon.
"""
from __future__ import annotations

import queue
import threading
import uuid
from datetime import datetime, timezone

import matlab_bridge as mb


_jobs = {}
_queues = {}
_worker_lock = threading.Lock()

_ALGO_DISPATCH = {
    "fseof":    mb.run_fseof,
    "optknock": mb.run_optknock,
    "optforce": mb.run_optforce,
    "oko":      mb.run_oko,
    "okoplus":  mb.run_oko_plus,
    "okoplus-build": mb.build_oko_intervals_from_homologs,
}


def _now_iso():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def create_job(model_id, biomass, target, algo, params):
    jid = uuid.uuid4().hex[:12]
    _jobs[jid] = {
        "status": "queued",
        "result": None,
        "error": None,
        "started_at": None,
        "ended_at": None,
        "track": "design",
        "algo": algo,
        "model_id": model_id,
        "target": target,
        "biomass": biomass,
        "params": dict(params or {}),
        "artifacts": {},
    }
    _queues[jid] = queue.Queue()
    threading.Thread(
        target=_run,
        args=(jid, model_id, biomass, target, algo, dict(params or {})),
        daemon=True,
    ).start()
    return jid


def cancel_job(jid):
    j = _jobs.get(jid)
    if not j:
        return False
    j["status"] = "cancelled"
    j["cancel_requested"] = True
    return True


def get_workspace_file(jid, name):
    j = _jobs.get(jid)
    if not j:
        return None
    return (j.get("artifacts") or {}).get(name)


def _emit(jid, ev):
    _queues[jid].put(ev)


def _bridge_args(algo, model_id, biomass, target, params):
    if algo == "fseof":
        return (model_id,
                int(params.get("iters", params.get("Iterations", 10))),
                float(params.get("coeff", params.get("Coefficient", 0.9))),
                biomass, target)
    if algo == "optknock":
        return (model_id, target, biomass, {
            "max_candidates":      int(params.get("max_candidates", 200)),
            "num_del":             int(params.get("num_del", 5)),
            "min_growth_fraction": float(params.get("min_growth_fraction", 0.1)),
            "vmax":                float(params.get("vmax", 1000)),
        })
    if algo == "optforce":
        return (model_id, target, biomass, {
            "k":              int(params.get("k", 2)),
            "nsets":          int(params.get("nsets", 1)),
            "max_candidates": int(params.get("max_candidates", 500)),
        })
    if algo == "oko":
        return (model_id, target, biomass, {
            "profile": params.get("profile", "auto"),
        })
    if algo == "okoplus":
        interval_path = params.get("interval_path", "")
        if not interval_path:
            raise mb.BridgeContractError(
                "err_param_invalid",
                "OKO+ requires an interval CSV path (rxn/uniprot/min/max).",
            )
        return (model_id, target, biomass, {
            "profile":       params.get("profile", "auto"),
            "interval_path": interval_path,
        })
    if algo == "okoplus-build":
        predictors = params.get("predictors") or ["UniKP"]
        max_homologs = int(params.get("max_homologs", params.get("maxHomologs", 100)))
        return (model_id, list(predictors), params.get("manager_path", ""), max_homologs)
    raise mb.BridgeContractError(
        "err_param_invalid",
        "unknown algorithm: " + str(algo),
    )


def _run(jid, model_id, biomass, target, algo, params):
    _jobs[jid]["status"] = "running"
    _jobs[jid]["started_at"] = _now_iso()
    _emit(jid, {"type": "log", "line": "queued->running (algo=" + str(algo) + ")"})

    if algo not in _ALGO_DISPATCH:
        _jobs[jid]["status"] = "error"
        _jobs[jid]["error"] = "unsupported algorithm: " + str(algo)
        _jobs[jid]["ended_at"] = _now_iso()
        _emit(jid, {"type": "error", "message": _jobs[jid]["error"]})
        return

    try:
        with _worker_lock:
            bridge_fn = _ALGO_DISPATCH[algo]
            try:
                args = _bridge_args(algo, model_id, biomass, target, params)
            except mb.BridgeContractError as pre:
                _jobs[jid]["status"] = "error"
                _jobs[jid]["error"] = pre.error_code + ": " + pre.error_message
                _jobs[jid]["ended_at"] = _now_iso()
                _emit(jid, {"type": "error", "message": _jobs[jid]["error"]})
                return
            out = bridge_fn(*args)
        # The homolog builder intentionally exposes the MATLAB envelope so
        # callers can inspect transport errors.  Jobs, however, have their own
        # status/error channel, so unwrap it here.
        if isinstance(out, dict) and "ok" in out:
            if not out.get("ok"):
                raise mb.BridgeContractError(
                    out.get("error_code", "err_internal"),
                    out.get("error_message", "OKO+ interval build failed"),
                )
            out = out.get("result") or {}
        _jobs[jid]["result"] = out
        if isinstance(out, dict):
            artifacts = {}
            for predictor, path in (out.get("predictor_csv_paths") or {}).items():
                artifacts[f"intervals_{str(predictor).lower()}.csv"] = path
            for predictor, path in (out.get("predictor_input_paths") or {}).items():
                artifacts[f"input_{str(predictor).lower()}.csv"] = path
            if out.get("candidate_csv_path"):
                artifacts["candidates.csv"] = out["candidate_csv_path"]
            if out.get("run_dir"):
                artifacts["run_dir"] = out["run_dir"]
            _jobs[jid]["artifacts"].update(artifacts)
        _jobs[jid]["status"] = "done"
        _jobs[jid]["ended_at"] = _now_iso()
        _emit(jid, {"type": "log", "line": "algo=" + str(algo) + " completed"})
        _emit(jid, {"type": "done"})
    except Exception as e:
        _jobs[jid]["status"] = "error"
        _jobs[jid]["error"] = str(e)
        _jobs[jid]["ended_at"] = _now_iso()
        _emit(jid, {"type": "error", "message": str(e)})


def get_job(jid):
    return _jobs.get(jid)


def get_queue(jid):
    return _queues.get(jid)


def get_result(jid):
    return _jobs.get(jid, {}).get("result")


def list_jobs(limit=3):
    if limit <= 0:
        return []
    items = []
    for jid, rec in _jobs.items():
        items.append({
            "jobId": jid,
            "status": rec.get("status", "?"),
            "track": rec.get("track"),
            "algo":  rec.get("algo"),
            "started_at": rec.get("started_at"),
            "ended_at": rec.get("ended_at"),
            "error": rec.get("error"),
        })
    items.sort(key=lambda x: (x.get("started_at") or "", x["jobId"]),
               reverse=True)
    return items[:limit]
