FROM rocker/r-ver@sha256:75020f8b3d011c81b98d05c47fc9234753b19ac64de7bf6ad2490adc1bf44eb9

ARG DANDELION_VERSION=0.1.0
ARG DANDELION_SHA256=8c1af5dae209f249c7a5f8700ac698c94a28d654ecb91ea038be9acc992ff69a
ENV DEBIAN_FRONTEND=noninteractive \
    R_REPOS=https://cloud.r-project.org \
    DANDELION_VERSION=${DANDELION_VERSION}

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 python3-pip ca-certificates curl \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --break-system-packages --no-cache-dir pyzipper==0.3.6 \
    && install2.r --error --skipinstalled jsonlite data.table igraph \
    && (curl -fsSL "https://cran.r-project.org/src/contrib/DANDELION_${DANDELION_VERSION}.tar.gz" -o /tmp/DANDELION.tar.gz \
        || curl -fsSL "https://cran.r-project.org/src/contrib/Archive/DANDELION/DANDELION_${DANDELION_VERSION}.tar.gz" -o /tmp/DANDELION.tar.gz) \
    && echo "${DANDELION_SHA256}  /tmp/DANDELION.tar.gz" | sha256sum -c - \
    && R CMD INSTALL /tmp/DANDELION.tar.gz \
    && rm -f /tmp/DANDELION.tar.gz \
    && Rscript -e "stopifnot(as.character(packageVersion('DANDELION')) == '0.1.0', !requireNamespace('qvalue', quietly = TRUE))"

RUN useradd --system --uid 10002 --create-home --shell /usr/sbin/nologin worker \
    && install -d -o worker -g worker /opt/nophigene

COPY --chown=worker:worker worker.py runner.R /opt/nophigene/

USER worker
WORKDIR /opt/nophigene
ENTRYPOINT ["python3", "/opt/nophigene/worker.py"]
