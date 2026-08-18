suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop(sprintf("%s is required", name))
  value
}

idat_prefix <- required_env("NOPHIGENE_IDAT_PREFIX")
output_dir <- required_env("NOPHIGENE_OUTPUT_DIR")
pipeline <- tolower(Sys.getenv("NOPHIGENE_METHYLATION_PIPELINE", unset = "sesame"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result <- tryCatch({
  if (pipeline == "sesame") {
    suppressPackageStartupMessages({library(sesame); library(sesameData)})
    sset <- openSesame(idat_prefix, prep = "QCDPB")
    beta <- getBetas(sset)
    detection <- tryCatch(pOOBAH(sset), error = function(e) rep(NA_real_, length(beta)))
    probes <- data.table(
      probe_id = names(beta), beta_value = as.numeric(beta), detection_p = as.numeric(detection),
      qc_pass = is.na(detection) | detection <= 0.01,
      qc_reason = ifelse(is.na(detection) | detection <= 0.01, "", "detection_p>0.01")
    )
    normalization <- "SeSAMe QCDPB"
  } else if (pipeline == "minfi_noob") {
    suppressPackageStartupMessages(library(minfi))
    targets <- data.frame(Basename = idat_prefix, stringsAsFactors = FALSE)
    rg <- read.metharray.exp(targets = targets, extended = TRUE)
    mset <- preprocessNoob(rg)
    beta_matrix <- getBeta(mset)
    detection_matrix <- detectionP(rg)
    probes <- data.table(
      probe_id = rownames(beta_matrix), beta_value = as.numeric(beta_matrix[, 1]),
      detection_p = as.numeric(detection_matrix[, 1]),
      qc_pass = as.numeric(detection_matrix[, 1]) <= 0.01,
      qc_reason = ifelse(as.numeric(detection_matrix[, 1]) <= 0.01, "", "detection_p>0.01")
    )
    normalization <- "minfi preprocessNoob"
  } else {
    stop("NOPHIGENE_METHYLATION_PIPELINE must be sesame or minfi_noob")
  }
  probes[, m_value := fifelse(beta_value > 0 & beta_value < 1, log2(beta_value / (1 - beta_value)), NA_real_)]
  probes[, normalization := normalization]
  fwrite(probes, file.path(output_dir, "measurements.csv"))
  list(
    ok = TRUE, pipeline = pipeline, normalization = normalization,
    probe_count = nrow(probes), qc_pass_count = sum(probes$qc_pass, na.rm = TRUE),
    detection_threshold = 0.01,
    policy = "Flagged probes are retained but excluded from primary statistics."
  )
}, error = function(error) {
  list(ok = FALSE, pipeline = pipeline, error = conditionMessage(error))
})

write_json(result, file.path(output_dir, "qc.json"), auto_unbox = TRUE, pretty = TRUE)
if (!isTRUE(result$ok)) quit(status = 1)
