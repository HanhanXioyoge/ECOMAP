# scripts/web/server/run.py
import os, sys, atexit, uvicorn
from pathlib import Path
import matlab_bridge as mb
from i18n_check import assert_parity, I18nParityError


# i18n dictionaries live at <repo>/scripts/web/matlab/i18n/{zh,en}.json.
# run.py lives at <repo>/scripts/web/server/run.py, so the i18n dir is three
# levels up from this file (server -> web -> scripts -> repo root), then back
# down into scripts/web/matlab/i18n.
_REPO_ROOT = Path(__file__).resolve().parents[3]
_I18N_DIR = _REPO_ROOT / "scripts" / "web" / "matlab" / "i18n"


def _resolve_i18n_paths():
    """Return ``(zh_path, en_path)`` for the parity check.

    Extracted so tests can ``monkeypatch.setattr`` it and point the server
    at a temp directory without touching the real dictionaries.
    """
    return _I18N_DIR / "zh.json", _I18N_DIR / "en.json"


def _write_pid_file():
    """Record our own PID so launchers can stop the whole tree reliably.

    The launcher cannot always learn the child PID itself (MATLAB ships a
    Java 8 JVM with no Process.pid(), and `pid` there resolves to the Control
    System Toolbox function), so the server is the source of truth.
    Path comes from MDP_PID_FILE, else <repo root>/.web.pid.
    """
    target = os.environ.get("MDP_PID_FILE")
    if target:
        pid_file = Path(target)
    else:
        # run.py lives at <root>/scripts/web/server/run.py
        pid_file = Path(__file__).resolve().parents[3] / ".web.pid"
    try:
        pid_file.write_text(str(os.getpid()), encoding="utf-8")
    except OSError:
        return  # not fatal: launcher falls back to port-based control
    atexit.register(lambda: pid_file.unlink(missing_ok=True))


def main(argv=None):
    """Server entry point.

    Performs the i18n parity check first so the server refuses to bind the
    port when zh.json and en.json have drifted apart. argv is accepted for
    future CLI parsing but currently unused.
    """
    zh_path, en_path = _resolve_i18n_paths()
    try:
        assert_parity(zh_path, en_path)
    except I18nParityError as exc:
        sys.stderr.write(
            f"[i18n] parity check failed: {exc}\n"
            f"[i18n] refusing to start; fix {zh_path} and {en_path} "
            f"so their key sets match.\n"
        )
        sys.exit(1)

    _write_pid_file()
    mb.get_engine()   # warm the MATLAB session before serving
    port = int(os.environ.get("MDP_PORT", "8000"))
    uvicorn.run("app:app", host="127.0.0.1", port=port, app_dir=os.path.dirname(__file__))


if __name__ == "__main__":
    main()