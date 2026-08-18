FROM bioconductor/bioconductor_docker:RELEASE_3_21

RUN R -e "BiocManager::install(c('sesame','sesameData','minfi'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('jsonlite','data.table'), repos='https://cloud.r-project.org')"

COPY scripts/methylation_preprocess.R /opt/nophigene/methylation_preprocess.R
ENTRYPOINT ["Rscript", "/opt/nophigene/methylation_preprocess.R"]
