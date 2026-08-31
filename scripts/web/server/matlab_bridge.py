"""
Python wrapper around the MATLAB Engine for the ECOMAP Web UI.

Every MATLAB response is checked against the envelope contract documented in
``scripts/web/matlab/bridge/CONTRACT.md`` before its result reaches a caller.
The canonical bridge failure codes are:

- ``err_init_fail``
- ``err_param_invalid``
- ``err_model_format``
- ``err_no_biomass``
- ``err_no_target``
- ``err_docker_missing``
- ``err_no_proteomics``
- ``err_gurobi_license``
- ``err_raven_notfound``
- ``err_sluice_data``
- ``err_kcat_merge``
- ``err_presto_data``
- ``err_oom``
- ``err_cancelled``

The public ``model_info``, ``algorithms``, and ``run_fseof`` helpers preserve
their existing result types. Contract violations and represented MATLAB
failures raise ``BridgeContractError`` with machine-readable error attributes.

The MATLAB engine is single-threaded and process-global, so this module lazily
starts one engine per Python process and serialises all calls with a lock.
RAVEN, COBRA, and other model toolboxes are expected on MATLAB's path. Optional
``MDP_RAVEN_DIR`` and ``MDP_COBRA_DIR`` environment variables are added
best-effort when present; no machine-specific paths are hard-coded.

UTF-8 source; LF line endings.
"""
from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from typing import Any

import matlab.engine


_BRIDGE_ERROR_CODES: frozenset[str] = frozenset(
    {
        "err_init_fail",
        "err_param_invalid",
        "err_model_format",
        "err_no_biomass",
        "err_no_target",
        "err_docker_missing",
        "err_no_proteomics",
        "err_gurobi_license",
        "err_raven_notfound",
        "err_sluice_data",
        "err_kcat_merge",
        "err_presto_data",
        "err_oom",
        "err_cancelled",
    }
)


class BridgeContractError(RuntimeError):
    """A malformed bridge response or a failure reported by MATLAB."""

    def __init__(self, error_code: str, error_message: str) -> None:
        self.error_code = error_code
        self.error_message = error_message
        super().__init__(f"{error_code}: {error_message}")

# ---------------------------------------------------------------------------
# Module state
# ---------------------------------------------------------------------------

_engine: "matlab.engine.MatlabEngine | None" = None
_engine_lock = threading.Lock()
_path_added = False

# Repo root = great-grandparent of this file
# (scripts/web/server/matlab_bridge.py
#  -> scripts/web/server -> scripts/web -> scripts -> repo).
# The MATLAB path needs to cover BOTH:
#   - scripts/             — shared core: loadModel, sniffModelType, algorithms/*
#   - scripts/web/matlab/  — web glue: bridge/, i18n/, registry/, result/
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_SRC_DIR = _REPO_ROOT / "scripts"
_WEB_MATLAB_DIR = _REPO_ROOT / "scripts" / "web" / "matlab"

# Optional toolbox directories, read from the environment. If set and present
# on disk, they are added to the MATLAB path. Missing/unset -> silent no-op.
# RAVEN / COBRA are otherwise expected to be on MATLAB's path already.
_OPTIONAL_TOOLBOX_ENV_VARS: tuple[str, ...] = ("MDP_RAVEN_DIR", "MDP_COBRA_DIR")


def _optional_toolbox_dirs() -> list[Path]:
    dirs: list[Path] = []
    for var in _OPTIONAL_TOOLBOX_ENV_VARS:
        val = os.environ.get(var)
        if val:
            dirs.append(Path(val))
    return dirs


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_envelope(env: Any, *, public_name: str) -> dict:
    """Validate and return one decoded MATLAB bridge envelope unchanged."""
    if not isinstance(env, dict):
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} returned {type(env).__name__}; expected a dict envelope",
        )

    required_keys = {"ok", "error_code", "error_message", "result"}
    missing_keys = required_keys.difference(env)
    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} response is missing required key(s): {missing}",
        )

    if not isinstance(env["ok"], bool):
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} response field 'ok' must be a bool",
        )
    if not isinstance(env["error_code"], str):
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} response field 'error_code' must be a str",
        )
    if not isinstance(env["error_message"], str):
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} response field 'error_message' must be a str",
        )
    if not env["ok"] and env["error_code"] not in _BRIDGE_ERROR_CODES:
        raise BridgeContractError(
            "bridge_contract_error",
            f"{public_name} returned unknown error code: {env['error_code']!r}",
        )

    return env


def _start_engine_locked() -> matlab.engine.MatlabEngine:
    """Start (if needed) and configure the engine. Caller MUST hold _engine_lock."""
    global _engine, _path_added
    if _engine is None:
        _engine = matlab.engine.start_matlab()
    if not _path_added:
        # Repository code first: shared core (scripts/) + web glue (scripts/web/matlab/).
        gp = _engine.genpath(str(_SRC_DIR))
        _engine.addpath(gp, nargout=0)
        gp = _engine.genpath(str(_WEB_MATLAB_DIR))
        _engine.addpath(gp, nargout=0)
        # Optional toolboxes from env vars, if installed (best-effort).
        for tb in _optional_toolbox_dirs():
            if tb.is_dir():
                try:
                    gp = _engine.genpath(str(tb))
                    _engine.addpath(gp, nargout=0)
                except Exception:  # noqa: BLE001
                    # Adding a path is never fatal; a missing/broken
                    # toolbox should not brick the engine.
                    pass
        _path_added = True
    return _engine


def _start_engine() -> matlab.engine.MatlabEngine:
    """Start a MATLAB engine, add scripts/ to its path, cache it globally."""
    with _engine_lock:
        return _start_engine_locked()


def _call_json(func_name: str, *args: Any) -> dict:
    """Call MATLAB, validate its decoded envelope, and raise on failure.

    Engine start and the MATLAB call are performed atomically under a single
    lock acquisition so no other thread can interleave between them.
    """
    with _engine_lock:
        eng = _start_engine_locked()
        response = getattr(eng, func_name)(*args, nargout=1)
    # MATLAB returns `matlab.char` (numpy-friendly str). json.loads accepts it.
    raw = json.loads(str(response))
    env = _validate_envelope(raw, public_name=func_name)
    if not env["ok"]:
        raise BridgeContractError(env["error_code"], env["error_message"])
    return env


def run_matlab(func_name: str, *args: Any) -> dict:
    """Call a MATLAB bridge function and return its decoded envelope unchanged.

    Unlike the typed wrappers below (``init_project``, ``load_model``,
    ...), this helper does **not** unwrap ``result`` and does **not**
    raise on ``ok=False``; the caller is responsible for inspecting the
    envelope fields directly. The envelope contract is still enforced:
    a malformed response (missing keys, wrong ``ok`` type, unknown error
    code on failure) raises ``BridgeContractError``.

    Engine start and the MATLAB call are performed atomically under a
    single lock acquisition so no other thread can interleave between
    them.

    Args:
        func_name: MATLAB function name (e.g. ``"mdpBuildOkoIntervalsFromHomologs"``).
        *args: positional arguments forwarded to the MATLAB function.

    Returns:
        The decoded envelope dict ``{ok, error_code, error_message, result}``.
    """
    with _engine_lock:
        eng = _start_engine_locked()
        response = getattr(eng, func_name)(*args, nargout=1)
    raw = json.loads(str(response))
    return _validate_envelope(raw, public_name=func_name)


def _unwrap(raw: Any) -> Any:
    """Validate a raw MATLAB envelope and return the success ``result`` payload.

    Used by the thin bridge wrappers (mdpInitProject, mdpLoadModel, ...) which
    call ``_engine.feval(...)`` directly. ``_call_json`` already returns the
    envelope; this helper is for wrappers that want to test-mock
    ``_engine.feval`` and validate the returned envelope themselves.
    """
    env = _validate_envelope(raw, public_name="_unwrap")
    if not env["ok"]:
        raise BridgeContractError(env["error_code"], env["error_message"])
    return env["result"]


def _feval(func_name: str, *args: Any, **kwargs: Any) -> Any:
    """Call ``func_name`` on the engine, lazily starting it on first use.

    The bridge wrappers used to read the module-global ``_engine`` directly,
    which silently relied on ``run.py`` having warmed it. Importing ``app:app``
    via any other entry point (test harness, uvicorn --reload, alternate
    launcher) would leave ``_engine is None`` and every feval would raise
    AttributeError. Routing through ``_start_engine()`` makes each call
    self-sufficient.

    Forwards kwargs (typically ``nargout=...``) to ``eng.feval``. We do NOT
    default nargout here because the per-call site already supplies it
    explicitly; mixing the two would risk a duplicate-kwarg TypeError.
    """
    eng = _start_engine()
    return eng.feval(func_name, *args, **kwargs)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_engine() -> matlab.engine.MatlabEngine:
    """Return the cached MATLAB engine, starting it on first call.

    Useful for callers (e.g. server/run.py) that want to warm the engine
    before serving HTTP traffic so the first request doesn't pay the
    MATLAB startup cost.
    """
    return _start_engine()


def model_info(path: str | os.PathLike[str]) -> dict:
    """Load a COBRA model from `path` and return a dict with metadata."""
    p = str(Path(path).resolve())
    envelope = _call_json("mdpModelInfo", p)
    return envelope["result"]


def algorithms(model_type: str) -> list:
    """Return the algorithm registry entries applicable to `model_type`.

    The MATLAB bridge ``mdpAlgorithms(modelType)`` already filters the
    registry server-side by checking ``modelType in entry.supports``.
    """
    envelope = _call_json("mdpAlgorithms", str(model_type))
    return envelope["result"]


def algorithms_by_track(track: str) -> list:
    """Return the algorithm registry entries whose ``track`` matches.

    Valid values are ``'recon'``, ``'calib'``, ``'analysis'``, ``'design'``,
    or ``'all'`` (no filtering). Implemented by the ``mdpAlgorithmsByTrack``
    MATLAB bridge on the server side; no client-side filtering is performed.
    """
    envelope = _call_json("mdpAlgorithmsByTrack", str(track))
    return envelope["result"]


def run_fseof(
    path: str | os.PathLike[str],
    biomass: str,
    target: str,
    iters: int,
    coeff: float,
) -> dict:
    """Run FSEOF and return {columns, rows, log} as a dict."""
    p = str(Path(path).resolve())
    # Positional args match scripts/web/matlab/bridge/mdpRunFseof.m:
    #   mdpRunFseof(path, biomassRxn, targetRxn, iterations, coefficient)
    envelope = _call_json("mdpRunFseof", p, str(biomass), str(target), int(iters), float(coeff))
    return envelope["result"]


# ---------------------------------------------------------------------------
# Reconstruction track (T6, T7, T8, T10, T11, T12, T13, T14)
# ---------------------------------------------------------------------------

def init_project(project_name: str, project_path: str) -> dict:
    """Initialise a new ECOMAP project; return {project_id, param_template_path, default_params}."""
    raw = _feval("mdpInitProject", project_name, project_path, nargout=1)
    return _unwrap(raw)


def load_parameter_manager(manager_path: str | os.PathLike[str]) -> dict:
    """Load and cache a project ParameterManager inside the MATLAB engine."""
    p = str(Path(manager_path).resolve())
    matlab_path = p.replace("'", "''")
    with _engine_lock:
        eng = _start_engine_locked()
        eng.eval(f"ParameterManager.getParams('{matlab_path}');", nargout=0)
    return {"managerPath": p}


def load_model(
    file_path: str,
    model_type: str = "Tradition",
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """Load a metabolic model; return {model_id, type, rxns, mets, genes, detected_organism}."""
    if manager_path:
        raw = _feval("mdpLoadModel", file_path, model_type, str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpLoadModel", file_path, model_type, nargout=1)
    return _unwrap(raw)


def convertec_model(
    model_id: str,
    topology: str,
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """S-matrix expansion; topology is one of basic/isozyme/integrated/all."""
    if manager_path:
        raw = _feval("mdpConvertecModel", model_id, topology, str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpConvertecModel", model_id, topology, nargout=1)
    return _unwrap(raw)


def annotate(
    ec_model_ids: list,
    annotation_stages: list,
    manager_path: str | os.PathLike[str] | None = None,
    options: dict | None = None,
) -> dict:
    """Annotate ecModels with protein complex / EC / metabolite info."""
    opts = options or {}
    if manager_path:
        raw = _feval("mdpAnnotate", list(ec_model_ids), list(annotation_stages), str(Path(manager_path).resolve()), opts, nargout=1)
    elif opts:
        raw = _feval("mdpAnnotate", list(ec_model_ids), list(annotation_stages), "", opts, nargout=1)
    else:
        raw = _feval("mdpAnnotate", list(ec_model_ids), list(annotation_stages), nargout=1)
    return _unwrap(raw)


def dl_predict(
    ec_model_id: str,
    models: list,
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """Generate DL kcat prediction inputs and execute predictions (Docker required)."""
    if manager_path:
        raw = _feval("mdpDlPredict", ec_model_id, list(models), str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpDlPredict", ec_model_id, list(models), nargout=1)
    return _unwrap(raw)


def kcat_compare(
    ec_model_id: str,
    dl_models: list,
    complex_names: list,
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """Compare DL-predicted kcat values against database kcat values."""
    if manager_path:
        raw = _feval("mdpKcatCompare", ec_model_id, list(dl_models), list(complex_names), str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpKcatCompare", ec_model_id, list(dl_models), list(complex_names), nargout=1)
    return _unwrap(raw)


def kcat_merge(
    ec_model_ids: list,
    use_custom_file: bool,
    options: dict | None = None,
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """Merge kcat from three sources (DB / DL / median) into all three ecModel types."""
    opts = options or {}
    if manager_path:
        raw = _feval("mdpKcatMerge", list(ec_model_ids), bool(use_custom_file), opts, str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpKcatMerge", list(ec_model_ids), bool(use_custom_file), opts, nargout=1)
    return _unwrap(raw)


def growth_predict(
    ec_model_id: str,
    c_source: str,
    bio_rxn: str,
    manager_path: str | os.PathLike[str] | None = None,
) -> dict:
    """Predict growth rate under three conditions; save ecModel to disk."""
    if manager_path:
        raw = _feval("mdpGrowthPredict", ec_model_id, c_source, bio_rxn, str(Path(manager_path).resolve()), nargout=1)
    else:
        raw = _feval("mdpGrowthPredict", ec_model_id, c_source, bio_rxn, nargout=1)
    return _unwrap(raw)


# ---------------------------------------------------------------------------
# Calibration track (T15, T16, T17, T18, T19, T20)
# ---------------------------------------------------------------------------

def apply_sluice(ec_model_id: str, ex_rxn_list: list) -> dict:
    """Apply the Sluice structure (prerequisite for all calibration methods)."""
    raw = _feval("mdpApplySluice", ec_model_id, list(ex_rxn_list), nargout=1)
    return _unwrap(raw)


def kcat_repo_init(ec_model_id: str) -> dict:
    """Initialise a fresh KcatRepo with the model's initial kcat set as 'Init' group."""
    raw = _feval("mdpKcatRepoInit", ec_model_id, nargout=1)
    return _unwrap(raw)


def sensitivity_tuning(ec_model_id: str, glc_ex: str,
                       target_growth: float, factor: float, multi: bool) -> dict:
    """Run single-condition (and optionally multi-condition) sensitivity analysis."""
    raw = _feval("mdpSensitivityTuning", ec_model_id, glc_ex,
                        float(target_growth), float(factor), bool(multi), nargout=1)
    return _unwrap(raw)


def gauks(ec_model_id: str, gem_model_id: str, bio_rxn: str) -> dict:
    """Run GAUKS calibration."""
    raw = _feval("mdpGauks", ec_model_id, gem_model_id, bio_rxn, nargout=1)
    return _unwrap(raw)


def bayesian(ec_model_id: str, scenarios: list, max_iter: int,
             proc: int, num_per_gen: int, reject_num: float,
             run_gauks_after: bool) -> dict:
    """Run Bayesian calibration across one or more scenarios."""
    raw = _feval("mdpBayesian", ec_model_id, scenarios,
                        int(max_iter), int(proc), int(num_per_gen),
                        float(reject_num), bool(run_gauks_after), nargout=1)
    return _unwrap(raw)


def presto(ec_model_id: str, data_files: dict) -> dict:
    """Run PRESTO calibration."""
    raw = _feval("mdpPresto", ec_model_id,
                        str(data_files["proteomics"]),
                        str(data_files["growth"]),
                        str(data_files["total_protein"]), nargout=1)
    return _unwrap(raw)


# ---------------------------------------------------------------------------
# Analysis track (T21, T22, T23a)
# ---------------------------------------------------------------------------

def ecfva(ec_model_id: str, target_rxn: str, fraction: float) -> dict:
    """Run ecFVA (FVA with enzyme constraint) for one target reaction."""
    raw = _feval("mdpEcFva", ec_model_id, target_rxn, float(fraction), nargout=1)
    return _unwrap(raw)


def knockout(ec_model_id: str, gene_list, c_source: str) -> dict:
    """Run single-gene knockout screen."""
    raw = _feval("mdpKnockout", ec_model_id, list(gene_list), c_source, nargout=1)
    return _unwrap(raw)


def protein_analysis(ec_model_id: str, group_by: str) -> dict:
    """Aggregate protein-usage statistics by a chosen group (e.g. subsystem)."""
    raw = _feval("mdpProteinAnalysis", ec_model_id, group_by, nargout=1)
    return _unwrap(raw)


# ---------------------------------------------------------------------------
# Design track (T23b, T23c)
# ---------------------------------------------------------------------------

def run_optknock(model_id: str, target: str, biomass: str, params: dict) -> dict:
    """Run OptKnock; return {candidates, summary}."""
    raw = _feval(
        "mdpRunOptknock",
        model_id, target, biomass,
        int(params.get("max_candidates", 200)),
        int(params.get("num_del", 5)),
        float(params.get("min_growth_fraction", 0.1)),
        float(params.get("vmax", 1000)),
        nargout=1,
    )
    return _unwrap(raw)


def run_optforce(model_id: str, target: str, biomass: str, params: dict) -> dict:
    """Run optForce; return {candidates, summary}."""
    raw = _feval(
        "mdpRunOptforce",
        model_id, target, biomass,
        int(params.get("k", 2)),
        int(params.get("nsets", 1)),
        int(params.get("max_candidates", 500)),
        nargout=1,
    )
    return _unwrap(raw)


def run_oko(model_id: str, target: str, biomass: str, params: dict) -> dict:
    """Run OKO; return {columns, rows, diagnostics, config, profile}.

    OKO minimises the number of kcat modifications required to overproduce a
    target metabolite, holding the enzyme abundance close to its wild-type
    allocation (Razaghi-Moghadam et al. 2024, Eqs 1-15).
    """
    raw = _feval(
        "mdpRunOko",
        model_id, biomass, target,
        str(params.get("profile", "auto")),
        nargout=1,
    )
    return _unwrap(raw)


def run_oko_plus(model_id: str, target: str, biomass: str, params: dict) -> dict:
    """Run OKO+; return {columns, rows, abundanceColumns, abundanceRows, ...}.

    OKO+ restricts kcat changes to ranges predicted by external DL models and
    only allows enzyme abundance changes at a high cost (paper Methods p.15,
    weight w=10). ``params['interval_path']`` must point to a CSV produced
    by ``buildOkoIntervals`` (rxn/uniprot/min/max columns).
    """
    interval_path = params.get("interval_path", "")
    if not interval_path:
        raise BridgeContractError(
            "err_param_invalid",
            "OKO+ requires an interval CSV path (rxn/uniprot/min/max).",
        )
    raw = _feval(
        "mdpRunOkoPlus",
        model_id, biomass, target,
        str(interval_path),
        str(params.get("profile", "auto")),
        nargout=1,
    )
    return _unwrap(raw)


def build_oko_intervals(predictions, predictor: str = "") -> dict:
    """Build OKO+ kcat intervals from one or more predictor tables.

    ``predictions`` is a list of ``{rows: [...], columns: [...]}`` payloads
    that match ``buildOkoIntervals``'s accepted schema (reaction/enzyme/
    organism/kcat). When ``predictor`` is non-empty only rows matching that
    predictor name are kept. Returns ``{rows, columns}`` in the same shape.
    """
    payload = [dict(p) for p in predictions]
    if predictor:
        raw = _feval(
            "mdpBuildOkoIntervals", payload, str(predictor), nargout=1,
        )
    else:
        raw = _feval(
            "mdpBuildOkoIntervals", payload, nargout=1,
        )
    return _unwrap(raw)


def build_oko_intervals_from_homologs(
    model_id: str,
    predictors: list,
    manager_path: str = "",
    max_homologs: int = 100,
) -> dict:
    """Build OKO+ kcat intervals from UniProt homolog retrieval + DL prediction.

    Wrapper around the MATLAB bridge ``mdpBuildOkoIntervalsFromHomologs``.
    Unlike the typed wrappers above, this function returns the envelope
    dict unchanged (no unwrap, no raise on ``ok=False``) because the
    OKO+ pipeline is an end-to-end orchestrator that the caller wants
    to inspect verbatim.

    Args:
        model_id: ecModel identifier (e.g., ``'eciML1515'``).
        predictors: list of predictor names; must be a subset of
            ``{'UniKP', 'CatPred'}``. ``DLKcat`` is **not** supported here
            because it has native multi-organism output already handled
            by ``build_oko_intervals`` (legacy).
        manager_path: optional ParameterManager.m path; pass ``""`` to
            use the currently loaded project.

    Returns:
        Envelope dict ``{ok, error_code, error_message, result}``
        where ``result`` (on success) has shape::

            {
                "predictor_csv_paths": {
                    "UniKP": str,    # absolute CSV path
                    "CatPred": str,  # absolute CSV path
                },
                "n_candidates_per_enzyme": [
                    {"rxn": str, "nHomologs": int},
                    ...
                ],
                "elapsed_seconds": float,
            }

        On validation failure (unsupported predictor) the function
        short-circuits with ``ok=False`` and ``error_code='err_param_invalid'``
        without calling MATLAB.
    """
    allowed = {"UniKP", "CatPred"}
    bad = set(predictors) - allowed
    if bad:
        return {
            "ok": False,
            "error_code": "err_param_invalid",
            "error_message": f"Unsupported predictors: {sorted(bad)}",
            "result": None,
        }

    return run_matlab(
        "mdpBuildOkoIntervalsFromHomologs",
        model_id,
        list(predictors),
        manager_path,
        float(max_homologs),
    )


# ---------------------------------------------------------------------------
# Teardown helper (mainly for tests / clean shutdown)
# ---------------------------------------------------------------------------

def shutdown() -> None:
    """Quit the cached MATLAB engine, if any. Idempotent."""
    global _engine, _path_added
    with _engine_lock:
        if _engine is not None:
            try:
                _engine.quit()
            except Exception:  # noqa: BLE001 - best-effort cleanup
                pass
            _engine = None
            _path_added = False
