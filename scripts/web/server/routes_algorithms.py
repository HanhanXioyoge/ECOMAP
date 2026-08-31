"""
Algorithms router — T31: GET /api/algorithms?track=...

Replaces the older ``GET /api/algorithms?modelType=...`` route which used
a different query parameter and returned the registry as a single flat
list. The new endpoint keeps the same path but:

- accepts an optional ``track`` query parameter (``recon`` / ``calib`` /
  ``analysis`` / ``design``);
- without ``track`` returns the full registry as a list;
- with ``track`` returns only the entries whose ``track`` field matches.

The MATLAB registry itself is fetched via ``matlab_bridge.algorithms``.
"""
from __future__ import annotations

from fastapi import APIRouter

import matlab_bridge as mb

algorithms_router = APIRouter()


@algorithms_router.get("/api/algorithms")
def algorithms_list(track: str = ""):
    """Return the algorithm registry; optionally filtered by ``track``."""
    registry = mb.algorithms("GEM")
    # Defensive: older MATLAB versions may return a dict keyed by track.
    if isinstance(registry, dict):
        if track and track in registry:
            return registry[track]
        return registry
    if not isinstance(registry, list):
        return []
    if track:
        return [a for a in registry
                if isinstance(a, dict) and a.get("track") == track]
    return registry
