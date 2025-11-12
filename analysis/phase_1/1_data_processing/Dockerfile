# ---- Base image ----
FROM python:3.11-slim-bullseye

# Avoid interactive apt, and keep layer lean
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# System deps: MAFFT + basics
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        mafft \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip && \
    pip install \
        numpy \
        pandas \
        scikit-learn \
        screed \
        openpyxl

COPY prepare-entry.py /prepare-entry.py
COPY gather-pairs.py /gather-pairs.py
COPY metrics.py /metrics.py