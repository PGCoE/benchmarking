# ---- Base image ----
FROM python:3.11-slim-bullseye AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_PROGRESS_BAR=off

RUN apt-get update && apt-get install -y --no-install-recommends \
      tini ca-certificates procps bash mafft \
    && rm -rf /var/lib/apt/lists/*

ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} app && useradd -m -u ${UID} -g ${GID} app

WORKDIR /app

# ---- Builder ----
FROM base AS builder
COPY pyproject.toml /app/
COPY src /app/src

RUN python -m pip install --upgrade pip setuptools wheel \
 && mkdir -p /wheels \
 && python -m pip wheel --wheel-dir /wheels /app

# ---- Runtime ----
FROM base AS runtime

COPY --from=builder /wheels /wheels
RUN python -m pip install --no-index --find-links=/wheels pgcoe_vbd \
 && rm -rf /wheels

# ✅ Copy your app source into runtime so the file path exists
COPY --from=builder /app/src /app/src

EXPOSE 8501

USER app
WORKDIR /workspace

ENTRYPOINT ["tini", "-g", "--"]
CMD ["streamlit", "run", "/app/src/pgcoe_vbd/app.py", \
     "--server.address=0.0.0.0", \
     "--server.port=8501", \
     "--server.enableCORS=false", \
     "--server.enableXsrfProtection=false"]