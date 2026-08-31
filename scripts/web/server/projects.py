"""Local JSON-backed project registry for the ECOMAP web UI."""

from __future__ import annotations

import json
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
PROJECTS_DIR = REPO_ROOT / "projects"
PROJECTS_FILE = Path(__file__).resolve().parent / "_uploads" / "projects" / "projects.json"
MODULE_DIRS = ("models", "reconstruction", "calibration", "analysis", "design")
RECONSTRUCTION_STEPS = (
    {
        "id": "initialize",
        "title": "Initialize project",
        "requires": (),
        "outputs": ("project.json", "*ParameterManagement.m"),
    },
    {
        "id": "loadGem",
        "title": "Load GEM",
        "requires": ("initialize",),
        "outputs": ("models/*.xml", "models/*.sbml", "models/*.json", "models/*.mat"),
    },
    {
        "id": "convert",
        "title": "Convert to ecGEMs",
        "requires": ("loadGem",),
        "outputs": ("models/*-integrated.mat", "models/*-basic.mat", "models/*-isozyme.mat"),
    },
    {
        "id": "annotate",
        "title": "Apply annotations",
        "requires": ("convert",),
        "outputs": ("reconstruction/uniprot.tsv", "reconstruction/ComplexPortal.json", "reconstruction/metInfo.tsv"),
    },
    {
        "id": "predictKcat",
        "title": "Predict Deep Learning kcat",
        "requires": ("annotate",),
        "outputs": (
            "reconstruction/kcatData/DLKcat.csv",
            "reconstruction/kcatData/UniKP.csv",
            "reconstruction/kcatData/CatPred.csv",
        ),
    },
    {
        "id": "compareKcat",
        "title": "Compare kcat predictions",
        "requires": ("predictKcat",),
        "outputs": ("analysis/AnalyzeKcatMatches/*match.mat", "analysis/AnalyzeKcatMatches/Benchmark_*.png"),
    },
    {
        "id": "mergeKcat",
        "title": "Merge kcat values",
        "requires": ("compareKcat",),
        "outputs": ("models/*-integrated.mat", "models/*merged*.mat", "models/*kcat*.mat"),
    },
    {
        "id": "growthSave",
        "title": "Growth validation and save models",
        "requires": ("mergeKcat",),
        "outputs": ("models/*-integrated.mat", "models/*-basic.mat", "models/*-isozyme.mat", "analysis/*growth*.csv", "analysis/*growth*.mat", "models/*calibrated*.mat"),
    },
)
PREVIEW_EXTENSIONS = {".png", ".jpg", ".jpeg", ".gif", ".webp"}
PARAMETER_MANAGER_FIELDS = {
    "InitialModel",
    "modeltype",
    "sigma",
    "Ptot",
    "f",
    "org_name",
    "uniprot.type",
    "uniprot.ID",
    "uniprot.geneIDfield",
    "uniprot.reviewed",
    "taxonomicID",
    "c_source",
    "bioRxn",
    "PRESTO.runParallel",
    "PRESTO.ncpu",
    "PRESTO.nIter",
    "PRESTO.epsilon",
    "PRESTO.lambda",
    "PRESTO.theta",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _stage_for_model_type(model_type: str) -> str:
    return "Ready for Calibration" if model_type == "ecGEM" else "Ready for Reconstruction"


def _display_organism(data: dict[str, Any], params: dict[str, Any]) -> str:
    return str(params.get("org_name") or data.get("org_name") or data.get("organism") or params.get("organism") or "")


def _matlab_identifier(value: str) -> str:
    name = re.sub(r"\W", "_", value.strip())
    if not name:
        name = "ECOMAP_Project"
    if not re.match(r"[A-Za-z]", name[0]):
        name = f"x{name}"
    return name


def _matlab_string(value: str) -> str:
    return str(value).replace("'", "''")


def _format_matlab_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return repr(value)
    return f"'{_matlab_string(str(value))}'"


def _flatten_params(params: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    for key, value in (params or {}).items():
        dotted = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            flattened.update(_flatten_params(value, dotted))
        else:
            flattened[dotted] = value
    return flattened


def _unflatten_params(params: dict[str, Any]) -> dict[str, Any]:
    nested: dict[str, Any] = {}
    for key, value in (params or {}).items():
        if "." in str(key):
            _assign_nested(nested, str(key), value)
        elif isinstance(value, dict):
            nested[str(key)] = _unflatten_params(value)
        else:
            nested[str(key)] = value
    return nested


def _project_folder_name(name: str) -> str:
    folder = re.sub(r'[<>:"/\\|?*]+', "_", name).strip().strip(".")
    if not folder:
        raise ValueError("project name cannot be used as a folder name")
    return folder


def _write_parameter_manager_template(project: dict[str, Any]) -> None:
    project_dir = Path(project["projectDir"])
    manager = str(project.get("parameterManager") or f"{project['projectId']}ParameterManagement.m")
    manager_path = project_dir / manager
    if manager_path.exists():
        return
    template_path = REPO_ROOT / "scripts" / "ParameterManagement" / "Template.m"
    text = template_path.read_text(encoding="utf-8")
    text = text.replace("KEY_Template", str(project["parameterManagerFunction"]))
    text = text.replace("KEY_PATH", _matlab_string(str(REPO_ROOT)))
    text = text.replace("KEY_NAME", _matlab_string(str(project["projectId"])))
    manager_path.write_text(text, encoding="utf-8")


def _write_parameter_manager_params(project: dict[str, Any], params: dict[str, Any]) -> None:
    project_dir = Path(project["projectDir"])
    manager = str(project.get("parameterManager") or f"{project['projectId']}ParameterManagement.m")
    manager_path = project_dir / manager
    if not manager_path.exists():
        _write_parameter_manager_template(project)
    text = manager_path.read_text(encoding="utf-8")
    values = {
        key: value
        for key, value in _flatten_params(params).items()
        if key in PARAMETER_MANAGER_FIELDS
    }
    for key, value in values.items():
        pattern = re.compile(
            rf"(?P<prefix>^\s*obj\.params\.{re.escape(key)}\s*=\s*)[^;\n]*(?P<suffix>;\s*(?:%.*)?$)",
            re.MULTILINE,
        )
        text, count = pattern.subn(
            lambda match, val=_format_matlab_value(value): f"{match.group('prefix')}{val}{match.group('suffix')}",
            text,
            count=1,
        )
        if count == 0:
            insert = f"    obj.params.{key} = {_format_matlab_value(value)};\n"
            text = re.sub(r"\nend\s*$", f"\n{insert}end", text, count=1)
    manager_path.write_text(text, encoding="utf-8")


def _read_registry() -> list[dict[str, Any]]:
    if not PROJECTS_FILE.exists():
        return []
    with open(PROJECTS_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    if isinstance(data, dict):
        return list(data.get("projects") or [])
    return list(data or [])


def _write_registry(projects: list[dict[str, Any]]) -> None:
    os.makedirs(PROJECTS_FILE.parent, exist_ok=True)
    tmp = PROJECTS_FILE.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump({"projects": projects}, handle, ensure_ascii=False, indent=2)
    os.replace(tmp, PROJECTS_FILE)


def _parse_literal(value: str) -> Any:
    value = value.strip()
    if not value:
        return ""
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1].replace("''", "'")
    low = value.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    try:
        if any(ch in value for ch in ".eE"):
            return float(value)
        return int(value)
    except ValueError:
        return value


def _assign_nested(target: dict[str, Any], dotted_key: str, value: Any) -> None:
    parts = [part for part in dotted_key.split(".") if part]
    if not parts:
        return
    cursor = target
    for part in parts[:-1]:
        next_value = cursor.get(part)
        if not isinstance(next_value, dict):
            next_value = {}
            cursor[part] = next_value
        cursor = next_value
    cursor[parts[-1]] = value


def _read_parameter_manager_params(project_dir: Path, data: dict[str, Any]) -> dict[str, Any]:
    manager = str(data.get("parameterManager") or "")
    if not manager:
        candidates = sorted(project_dir.glob("*ParameterManagement.m"))
        if not candidates:
            return {}
        manager_path = candidates[0]
    else:
        manager_path = project_dir / manager
    if not manager_path.exists():
        return {}
    try:
        text = manager_path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = manager_path.read_text(encoding="utf-8", errors="ignore")
    params: dict[str, Any] = {}
    pattern = re.compile(r"obj\.params\.([A-Za-z]\w*(?:\.[A-Za-z]\w*)*)\s*=\s*([^;\n]+)\s*;")
    for match in pattern.finditer(text):
        key, raw_value = match.groups()
        if "fullfile(" in raw_value or "mfilename" in raw_value:
            continue
        _assign_nested(params, key, _parse_literal(raw_value))
    return params


def _normalise_project(data: dict[str, Any], project_dir: Path | None = None) -> dict[str, Any]:
    if not isinstance(data, dict):
        return {}
    project_id = str(
        data.get("projectId")
        or data.get("project_id")
        or data.get("projectName")
        or data.get("name")
        or (project_dir.name if project_dir else uuid.uuid4().hex[:12])
    )
    params = data.get("params") if isinstance(data.get("params"), dict) else {}
    if project_dir is not None:
        params = {**_read_parameter_manager_params(project_dir, data), **params}
    project = dict(data)
    project["projectId"] = project_id
    project["name"] = str(data.get("name") or data.get("projectName") or project_id)
    project["projectName"] = str(data.get("projectName") or project["name"])
    if project_dir is not None:
        project["projectDir"] = str(project_dir)
    elif data.get("projectDir"):
        project["projectDir"] = str(data["projectDir"])
    project["modelId"] = str(data.get("modelId") or data.get("InitialModel") or params.get("InitialModel") or "")
    input_model_type = str(data.get("modelType") or data.get("inputModelType") or "")
    if not input_model_type:
        input_model_type = "ecGEM" if str(params.get("modeltype") or "").lower() in {"ecomap", "smoment", "ecmpy", "gecko"} else "GEM"
    project["modelType"] = input_model_type
    project["organism"] = _display_organism(data, params)
    project["params"] = params
    project["stage"] = str(data.get("stage") or _stage_for_model_type(project["modelType"] or "GEM"))
    project["status"] = str(data.get("status") or "ready")
    project["createdAt"] = str(data.get("createdAt") or data.get("updatedAt") or "")
    project["updatedAt"] = str(data.get("updatedAt") or data.get("createdAt") or "")
    return project


def _read_local_projects() -> list[dict[str, Any]]:
    if not PROJECTS_DIR.exists():
        return []
    projects: list[dict[str, Any]] = []
    for project_dir in PROJECTS_DIR.iterdir():
        if not project_dir.is_dir():
            continue
        meta_path = project_dir / "project.json"
        if not meta_path.exists():
            continue
        try:
            with open(meta_path, encoding="utf-8") as handle:
                project = _normalise_project(json.load(handle), project_dir)
        except (OSError, json.JSONDecodeError):
            continue
        if project:
            projects.append(project)
    return projects


def _read_all() -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for project in _read_registry():
        normalised = _normalise_project(project)
        merged[normalised["projectId"]] = normalised
    for project in _read_local_projects():
        merged[project["projectId"]] = project
    return list(merged.values())


def _write_project_json(project: dict[str, Any]) -> None:
    project_dir = Path(project["projectDir"])
    os.makedirs(project_dir, exist_ok=True)
    for dirname in MODULE_DIRS:
        os.makedirs(project_dir / dirname, exist_ok=True)
    payload = dict(project)
    payload["directories"] = {
        dirname: str(project_dir / dirname)
        for dirname in MODULE_DIRS
    }
    tmp = project_dir / "project.tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(tmp, project_dir / "project.json")


def _require_text(payload: dict[str, Any], key: str) -> str:
    value = str(payload.get(key) or "").strip()
    if not value:
        raise ValueError(f"missing field: {key}")
    return value


def create_project(payload: dict[str, Any]) -> dict[str, Any]:
    name = _require_text(payload, "name")
    model_type = _require_text(payload, "modelType")
    if model_type not in {"GEM", "ecGEM"}:
        raise ValueError("modelType must be GEM or ecGEM")

    now = _now()
    folder_name = _project_folder_name(name)
    project_dir = PROJECTS_DIR / folder_name
    if (project_dir / "project.json").exists():
        raise ValueError("project already exists")
    manager_function = f"{_matlab_identifier(folder_name)}ParameterManagement"
    incoming_params = payload.get("params") if isinstance(payload.get("params"), dict) else {}
    params = _unflatten_params({
        **_flatten_params({
            "projectName": folder_name,
            "InitialModel": "",
            "modeltype": "Tradition",
        }),
        **_flatten_params(incoming_params),
    })
    project = {
        "projectId": folder_name,
        "name": name,
        "projectName": folder_name,
        "projectDir": str(project_dir),
        "modelId": "",
        "modelType": model_type,
        "organism": _display_organism(payload, params),
        "params": params,
        "stage": payload.get("stage") or _stage_for_model_type(model_type),
        "status": "ready",
        "createdAt": now,
        "updatedAt": now,
        "parameterManager": f"{manager_function}.m",
        "parameterManagerFunction": manager_function,
        "runs": [],
        "files": [],
    }
    _write_project_json(project)
    _write_parameter_manager_template(project)
    projects = _read_registry()
    projects.append(project)
    _write_registry(projects)
    return dict(project)


def list_projects(limit: int | None = None) -> list[dict[str, Any]]:
    projects = sorted(_read_all(), key=lambda item: item.get("updatedAt", ""), reverse=True)
    if limit is None:
        return projects
    return projects[: max(0, int(limit))]


def get_project(project_id: str) -> dict[str, Any] | None:
    for project in _read_all():
        if project.get("projectId") == project_id:
            return dict(project)
    return None


def _safe_project_dir(project_id: str) -> Path:
    project = get_project(project_id)
    if not project or not project.get("projectDir"):
        raise KeyError(project_id)
    project_dir = Path(project["projectDir"]).resolve()
    projects_root = PROJECTS_DIR.resolve()
    try:
        project_dir.relative_to(projects_root)
    except ValueError as exc:
        raise KeyError(project_id) from exc
    return project_dir


def _file_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in PREVIEW_EXTENSIONS:
        return "image"
    if suffix in {".mat"}:
        return "model"
    if suffix in {".xml", ".sbml", ".json", ".yml", ".yaml"}:
        return "data"
    if suffix in {".csv", ".tsv", ".txt"}:
        return "table"
    return "file"


def _scan_files(project_dir: Path) -> list[dict[str, Any]]:
    files: list[dict[str, Any]] = []
    for dirname in MODULE_DIRS:
        base = project_dir / dirname
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(project_dir).as_posix()
            stat = path.stat()
            files.append({
                "name": path.name,
                "relativePath": rel,
                "folder": rel.split("/", 1)[0],
                "kind": _file_kind(path),
                "size": stat.st_size,
                "updatedAt": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                "previewable": path.suffix.lower() in PREVIEW_EXTENSIONS,
            })
    return files


def get_project_files(project_id: str) -> dict[str, Any]:
    project_dir = _safe_project_dir(project_id)
    files = _scan_files(project_dir)
    groups: dict[str, list[dict[str, Any]]] = {}
    for item in files:
        groups.setdefault(item["folder"], []).append(item)
    return {
        "projectId": project_id,
        "projectDir": str(project_dir),
        "files": files,
        "groups": groups,
    }


def resolve_project_file(project_id: str, relative_path: str) -> Path:
    project_dir = _safe_project_dir(project_id)
    target = (project_dir / relative_path).resolve()
    try:
        target.relative_to(project_dir)
    except ValueError as exc:
        raise KeyError(relative_path) from exc
    if not target.is_file():
        raise KeyError(relative_path)
    return target


def project_model_upload_path(project_id: str, filename: str) -> Path:
    project_dir = _safe_project_dir(project_id)
    safe_name = os.path.basename(filename or "model.mat")
    if not safe_name:
        safe_name = "model.mat"
    models_dir = (project_dir / "models").resolve()
    os.makedirs(models_dir, exist_ok=True)
    target = (models_dir / safe_name).resolve()
    try:
        target.relative_to(models_dir)
    except ValueError as exc:
        raise KeyError(filename) from exc
    return target


def parameter_manager_path(project_id: str) -> Path:
    project = get_project(project_id)
    if not project:
        raise KeyError(project_id)
    project_dir = _safe_project_dir(project_id)
    manager = str(project.get("parameterManager") or "")
    candidates = [project_dir / manager] if manager else []
    candidates.extend(sorted(project_dir.glob("*ParameterManagement.m")))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise KeyError(project_id)


def _matches_any(project_dir: Path, patterns: tuple[str, ...]) -> list[str]:
    matched: list[str] = []
    for pattern in patterns:
        for path in project_dir.glob(pattern):
            if path.is_file():
                matched.append(path.relative_to(project_dir).as_posix())
    return sorted(set(matched))


def get_reconstruction_state(project_id: str) -> dict[str, Any]:
    project = get_project(project_id)
    if not project:
        raise KeyError(project_id)
    project_dir = _safe_project_dir(project_id)
    completed: set[str] = set()
    steps: list[dict[str, Any]] = []
    for spec in RECONSTRUCTION_STEPS:
        outputs = _matches_any(project_dir, tuple(spec["outputs"]))
        if spec["id"] == "initialize" and not str((project.get("params") or {}).get("InitialModel") or "").strip():
            outputs = []
        required = tuple(spec["requires"])
        dependencies_done = all(step_id in completed for step_id in required)
        if outputs and dependencies_done:
            status = "completed"
            completed.add(spec["id"])
        elif dependencies_done:
            status = "ready"
        else:
            status = "locked"
        steps.append({
            "id": spec["id"],
            "title": spec["title"],
            "status": status,
            "requires": list(required),
            "outputs": outputs,
        })
    return {
        "projectId": project_id,
        "project": project,
        "steps": steps,
        "files": get_project_files(project_id)["files"],
    }


def update_project_params(project_id: str, params: dict[str, Any]) -> dict[str, Any]:
    projects = _read_registry()
    target = get_project(project_id)
    if not target:
        raise KeyError(project_id)
    merged_flat = {
        **_flatten_params(target.get("params") or {}),
        **_flatten_params(params or {}),
    }
    merged = _unflatten_params(merged_flat)
    target["params"] = merged
    target["organism"] = _display_organism(target, merged)
    target["modelId"] = merged.get("InitialModel", target.get("modelId", ""))
    target["stage"] = target.get("stage") or _stage_for_model_type(target.get("modelType", "GEM"))
    target["updatedAt"] = _now()
    if target.get("projectDir"):
        _write_parameter_manager_params(target, merged)
        _write_project_json(target)

    found_registry = False
    for idx, project in enumerate(projects):
        if project.get("projectId") != project_id:
            continue
        projects[idx] = target
        found_registry = True
        break
    if not found_registry:
        projects.append(target)
    _write_registry(projects)
    return dict(target)
