# ---- Base image ----
FROM python:3.11-slim-bullseye

# Avoid interactive apt, and keep layer lean
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    # Streamlit defaults
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_HEADLESS=true \
    # Matplotlib cache path (prevents permission issues)
    MPLCONFIGDIR=/tmp/matplotlib

# System deps: MAFFT + basics
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        mafft \
        ca-certificates \
        tini \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user & workdir
RUN useradd -m -u 1000 appuser
WORKDIR /app

RUN pip install --upgrade pip && \
    pip install \
        numpy \
        pandas \
        scikit-learn \
        plotly \
        streamlit \
        screed \
        openpyxl

COPY prepare-entry.py /app/prepare-entry.py
COPY gather-pairs.py /app/gather-pairs.py
COPY metrics.py /app/metrics.py
COPY dashboard.py /app/dashboard.py

# Switch to non-root
USER appuser

# Streamlit port
EXPOSE 8501

# Use tini as PID 1 for clean signal handling
ENTRYPOINT ["/usr/bin/tini", "--"]

# Default to running the Streamlit app
CMD ["run_dashboard"]
