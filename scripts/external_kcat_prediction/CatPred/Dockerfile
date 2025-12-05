# =============================
# Base: CUDA 12.8 + cuDNN (for RTX 50 series / Blackwell)
# =============================
FROM nvidia/cuda:12.8.0-cudnn-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

# 1) System basics
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bash tini git \
      libxrender1 libxext6 libsm6 libx11-6 libxcb1 libx11-xcb1 \
      libcairo2 libfreetype6 libfontconfig1 \
 && rm -rf /var/lib/apt/lists/*

# 2) Miniforge
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN curl -fsSL -o /tmp/miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
 && bash /tmp/miniforge.sh -b -p "$CONDA_DIR" \
 && rm -f /tmp/miniforge.sh

# 3) Conda config
RUN conda config --set solver libmamba \
 && conda config --set channel_priority strict

WORKDIR /app

COPY environment.yml /tmp/environment.yml

# 4) Create env from the filtered yml
RUN conda env create -f /tmp/environment.yml \
 && conda clean -afy
#    把该 env 放到 PATH
ENV CONDA_DEFAULT_ENV=catpred
ENV PATH=$CONDA_DIR/envs/catpred/bin:$PATH
SHELL ["/bin/bash", "-lc"]

RUN python -m pip install --upgrade pip \
 && pip install --index-url https://download.pytorch.org/whl/cu128 torch==2.8.0 \
 && python - <<'PY'
import torch
print(">>> torch after pip:", torch.__version__, "cuda:", torch.version.cuda, "avail:", torch.cuda.is_available())
PY

RUN pip install --no-deps \
      fair-esm==2.0.0 \
      rotary-embedding-torch==0.6.5

RUN pip install \
      ipdb \
      progres \
      rdkit-pypi \
      transformers \
      descriptastorus \
      sentencepiece \
      seaborn \
      typed-argument-parser>=1.6.1

ENV TORCH_HOME=/root/.cache/torch
RUN mkdir -p $TORCH_HOME/hub/checkpoints
COPY esm2_t33_650M_UR50D.pt \
     $TORCH_HOME/hub/checkpoints/esm2_t33_650M_UR50D.pt

COPY . /app

RUN python - <<'PY'
import torch
print("FINAL:", torch.__version__, torch.version.cuda, torch.cuda.is_available())
if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
PY