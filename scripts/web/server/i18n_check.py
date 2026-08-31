"""
Lightweight i18n parity check.

The ECOMAP web frontend renders its UI from flat ``scripts/web/matlab/i18n/{zh,en}.json``
dictionaries. Both files must keep the same set of keys, otherwise some
text will render as a raw i18n key on one side or the other.

A missing key does NOT stop the application from running — the new design
(B10) explicitly chooses the ``inherit + extend`` policy (沿用 + 补), so an
entry may exist in one file but not the other for a transitional period.
``assert_parity`` is meant to be invoked in CI / pre-merge to catch drift.

UTF-8 source; LF line endings.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable


class I18nParityError(RuntimeError):
    """Raised when two i18n dictionaries do not share the same key set."""


def _load_dict(path: Path) -> dict:
    """Load a JSON file and return it as a dict. Empty file -> empty dict."""
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise I18nParityError(
            f"{path}: top-level JSON must be an object/dict, not {type(data).__name__}"
        )
    return data


def _diff_keys(a: Iterable[str], b: Iterable[str]) -> tuple[set[str], set[str]]:
    """Return ``(missing_in_a, missing_in_b)`` for two key iterables."""
    sa, sb = set(a), set(b)
    return sb - sa, sa - sb


def assert_parity(zh_path: str | Path, en_path: str | Path) -> None:
    """Assert that ``zh_path`` and ``en_path`` have the same set of keys.

    Raises :class:`I18nParityError` (a ``RuntimeError`` subclass) listing
    exactly which keys are present on one side but missing on the other.
    """
    zh_p, en_p = Path(zh_path), Path(en_path)
    zh = _load_dict(zh_p)
    en = _load_dict(en_p)

    missing_in_zh, missing_in_en = _diff_keys(zh.keys(), en.keys())
    if not missing_in_zh and not missing_in_en:
        return

    parts: list[str] = []
    if missing_in_zh:
        parts.append(f"keys missing in {zh_p}: {sorted(missing_in_zh)}")
    if missing_in_en:
        parts.append(f"keys missing in {en_p}: {sorted(missing_in_en)}")
    raise I18nParityError("; ".join(parts))


def get_all_keys(zh_path: str | Path, en_path: str | Path) -> list[str]:
    """Return the sorted union of keys present in either file."""
    zh_p, en_p = Path(zh_path), Path(en_path)
    zh = _load_dict(zh_p)
    en = _load_dict(en_p)
    return sorted(set(zh) | set(en))


__all__ = ["I18nParityError", "assert_parity", "get_all_keys"]
