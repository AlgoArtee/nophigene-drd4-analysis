FROM python:3.10-slim AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /build

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev && \
    rm -rf /var/lib/apt/lists/*

COPY requirements-app.txt .

RUN python -m pip install --upgrade pip && \
    python -m pip wheel --wheel-dir /wheels -r requirements-app.txt


FROM python:3.10-slim

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

ARG USERNAME=appuser
ARG USER_UID=1000
ARG USER_GID=1000

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 && \
    rm -rf /var/lib/apt/lists/*

RUN groupadd --gid "${USER_GID}" "${USERNAME}" && \
    useradd --uid "${USER_UID}" --gid "${USER_GID}" -m "${USERNAME}"

WORKDIR /home/${USERNAME}/app

COPY --from=builder /wheels /wheels
COPY requirements-app.txt .

RUN python -m pip install --upgrade pip && \
    python -m pip install --no-index --find-links=/wheels -r requirements-app.txt

COPY --chown=${USERNAME}:${USERNAME} src ./src

RUN mkdir -p data results && \
    chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/app

USER ${USERNAME}

EXPOSE 8000

ENTRYPOINT ["python", "src/app.py"]
CMD ["web", "--host", "0.0.0.0", "--port", "8000"]
