FROM python:3.12-slim-bookworm

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    NOPHIGENE_MODEL_WORKER_IMAGE=nophigene/alphagenome-runner:0.8.0

LABEL org.opencontainers.image.title="NophiGene AlphaGenome isolated worker" \
      org.opencontainers.image.version="0.8.0" \
      org.opencontainers.image.revision="71a6beb8c30832f121309a81c2530efa5af7986a"

RUN groupadd --system --gid 10003 modelworker \
    && useradd --system --uid 10003 --gid modelworker --create-home --shell /usr/sbin/nologin modelworker \
    && install -d -o modelworker -g modelworker /opt/nophigene /work /credentials

COPY requirements-sdk.txt /tmp/requirements-sdk.txt

# Install the reviewed SDK's declared runtime dependencies, then install the
# official wheel with PyPI's published SHA-256 and without dependency changes.
RUN python -m pip install --no-cache-dir \
      absl-py anndata fsspec 'grpcio>=1.67.1' immutabledict jaxtyping matplotlib \
      ml_dtypes numpy pandas 'protobuf>=5.28.3' pyarrow pyfaidx scipy seaborn \
      tqdm typeguard typing_extensions zstandard pyzipper==0.3.6 \
    && python -m pip install --no-cache-dir --no-deps --require-hashes -r /tmp/requirements-sdk.txt \
    && python -c "from importlib.metadata import version; assert version('alphagenome') == '0.8.0'"

COPY --chown=modelworker:modelworker prepare_request.py worker.py /opt/nophigene/

USER modelworker
WORKDIR /opt/nophigene
ENTRYPOINT ["python", "/opt/nophigene/worker.py"]
