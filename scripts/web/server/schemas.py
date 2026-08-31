# server/schemas.py
"""Pydantic v2 request/response models for the MDP Web UI backend."""
from pydantic import BaseModel


class ModelResp(BaseModel):
    modelId: str
    modelType: str
    nRxns: int
    nMets: int
    nGenes: int
    rxns: list[str]
