FROM nvidia/cuda:13.0.1-cudnn-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bash tini \
 && rm -rf /var/lib/apt/lists/*

RUN curl -L -o /tmp/miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
 && bash /tmp/miniforge.sh -b -p /opt/conda \
 && rm -f /tmp/miniforge.sh

RUN /opt/conda/bin/conda config --set channel_priority strict \
 && /opt/conda/bin/conda create -y -n py39 python=3.9 \
 && /opt/conda/bin/conda clean -afy
ENV PATH=/opt/conda/envs/py39/bin:/opt/conda/bin:$PATH

WORKDIR /app

COPY requirements.txt /app/requirements.txt

RUN python -m pip install --upgrade pip \
 && python -m pip install torch torchvision --index-url https://download.pytorch.org/whl/cu129 \
 && python -m pip install -r /app/requirements.txt

COPY DLKcat /app/DLKcat

COPY UniKP /app/UniKP

CMD ["python", "-c", "import torch, torchvision; print('torch=', torch.__version__, 'cuda=', torch.cuda.is_available())"]