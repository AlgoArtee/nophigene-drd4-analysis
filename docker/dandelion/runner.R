#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(DANDELION)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: runner.R REQUEST_JSON RESULT_JSON")
request_path <- args[[1L]]
result_path <- args[[2L]]
request <- jsonlite::fromJSON(request_path, simplifyVector = FALSE)

if (as.character(utils::packageVersion("DANDELION")) != "0.1.0") {
  stop("Allowlisted DANDELION package version 0.1.0 is required.")
}
if (requireNamespace("qvalue", quietly = TRUE)) {
  stop("qvalue must not be installed: this runner requires deterministic Benjamini-Hochberg adjustment.")
}

job_dir <- request$job_dir
cancel_file <- file.path(job_dir, "cancel.requested")
write_progress <- function(percent, stage) {
  jsonlite::write_json(
    list(progress_percent = as.integer(percent), stage = stage),
    file.path(job_dir, "progress.json"),
    auto_unbox = TRUE
  )
}
check_cancel <- function() {
  if (file.exists(cancel_file)) stop("Cancellation requested.")
}

file_spec <- function(role) {
  item <- request$dataset$files[[role]]
  if (is.null(item)) stop(sprintf("Missing file role: %s", role))
  item
}
mapping_name <- function(spec, canonical, fallback) {
  value <- spec$mapping[[canonical]]
  if (is.null(value) || !nzchar(as.character(value))) fallback else as.character(value)
}
read_tabular <- function(spec, select = NULL) {
  if (spec$format == "rds") return(readRDS(spec$runtime_path))
  separator <- if (spec$format == "csv") "," else "\t"
  data.table::fread(spec$runtime_path, sep = separator, select = select, data.table = FALSE, check.names = FALSE)
}
rename_required <- function(data, spec, definitions) {
  for (canonical in names(definitions)) {
    source <- mapping_name(spec, canonical, definitions[[canonical]])
    if (!(source %in% colnames(data))) stop(sprintf("%s is missing mapped column '%s'.", spec$role, source))
    colnames(data)[colnames(data) == source] <- canonical
  }
  data
}
validate_p <- function(values, label) {
  numbers <- suppressWarnings(as.numeric(values))
  invalid <- !is.na(numbers) & (!is.finite(numbers) | numbers < 0 | numbers > 1)
  if (any(invalid)) stop(sprintf("%s contains p-values outside [0, 1].", label))
  numbers
}

write_progress(2, "loading_reference_inputs")
association_spec <- file_spec("gene_association")
association <- read_tabular(association_spec)
association <- rename_required(association, association_spec, c(gene = "gene", p_value = "p_value"))
if (anyDuplicated(association$gene)) stop("gene_association contains duplicate gene identifiers.")
p_wes <- validate_p(association$p_value, "gene_association")
names(p_wes) <- as.character(association$gene)

annotation_spec <- file_spec("gene_annotation")
annotation <- read_tabular(annotation_spec)
annotation <- rename_required(
  annotation,
  annotation_spec,
  c(gene_name = "gene_name", type = "type", Chromosome = "Chromosome", start = "start", end = "end")
)
annotation$gene_name <- as.character(annotation$gene_name)
annotation$Chromosome <- as.character(annotation$Chromosome)
annotation$Chromosome <- ifelse(grepl("^chr", annotation$Chromosome), annotation$Chromosome, paste0("chr", annotation$Chromosome))
annotation$start <- as.numeric(annotation$start)
annotation$end <- as.numeric(annotation$end)

snp_ref <- NULL
if (request$dataset$exposure_type == "SNP") {
  snp_spec <- file_spec("snp_reference")
  snp_ref <- read_tabular(snp_spec)
  snp_ref <- rename_required(snp_ref, snp_spec, c(SNP = "SNP", SNPPos = "SNPPos", SNPChr = "SNPChr"))
}

trans_spec <- file_spec("trans_matrix")
if (trans_spec$format == "rds") {
  trans_all <- readRDS(trans_spec$runtime_path)
  if (is.data.frame(trans_all)) {
    row_column <- mapping_name(trans_spec, "row_id", colnames(trans_all)[1L])
    rownames(trans_all) <- as.character(trans_all[[row_column]])
    trans_all[[row_column]] <- NULL
    trans_all <- as.matrix(trans_all)
  }
  if (!is.matrix(trans_all)) stop("RDS trans_matrix must contain a numeric matrix or data frame.")
  exposure_names <- colnames(trans_all)
  row_id_column <- NULL
} else {
  separator <- if (trans_spec$format == "csv") "," else "\t"
  header <- data.table::fread(trans_spec$runtime_path, sep = separator, nrows = 0L, data.table = FALSE, check.names = FALSE)
  row_id_column <- mapping_name(trans_spec, "row_id", colnames(header)[1L])
  if (!(row_id_column %in% colnames(header))) stop("trans_matrix row_id column is missing.")
  exposure_names <- setdiff(colnames(header), row_id_column)
  trans_all <- NULL
}

if (!is.null(request$dataset$files$exposure_list)) {
  exposure_spec <- file_spec("exposure_list")
  exposure_data <- read_tabular(exposure_spec)
  exposure_column <- mapping_name(exposure_spec, "exposure", colnames(exposure_data)[1L])
  exposure_names <- intersect(exposure_names, as.character(exposure_data[[exposure_column]]))
}
if (length(exposure_names) == 0L) stop("No exposure columns are eligible for analysis.")

chunk_size <- as.integer(request$parameters$chunk_size)
chunks <- split(exposure_names, ceiling(seq_along(exposure_names) / chunk_size))
target_fdr <- as.numeric(request$parameters$target_fdr)
cis_window <- as.numeric(request$parameters$cis_window_bp)
gene_threshold <- as.numeric(request$parameters$gene_association_threshold)
all_records <- list()
record_index <- 0L
tested_pairs <- 0L
artifact_paths <- character()

snp_map <- NULL
if (!is.null(request$dataset$files$snp_gene_map)) {
  map_spec <- file_spec("snp_gene_map")
  snp_map <- read_tabular(map_spec)
  snp_map <- rename_required(snp_map, map_spec, c(SNP = "SNP", GeneSymbol = "GeneSymbol"))
}

for (chunk_index in seq_along(chunks)) {
  check_cancel()
  selected <- chunks[[chunk_index]]
  write_progress(5 + floor(88 * (chunk_index - 1L) / length(chunks)), sprintf("dandelion_chunk_%s_of_%s", chunk_index, length(chunks)))
  if (is.null(trans_all)) {
    frame <- read_tabular(trans_spec, select = c(row_id_column, selected))
    rownames(frame) <- as.character(frame[[row_id_column]])
    frame[[row_id_column]] <- NULL
    p_trans <- as.matrix(frame)
  } else {
    p_trans <- trans_all[, selected, drop = FALSE]
  }
  storage.mode(p_trans) <- "numeric"
  invalid <- !is.na(p_trans) & (!is.finite(p_trans) | p_trans < 0 | p_trans > 1)
  if (any(invalid)) stop("trans_matrix contains p-values outside [0, 1].")
  result <- DANDELION::med_gene(
    p.trans = p_trans,
    p.wes = p_wes,
    ref.table = annotation,
    gene1.list = selected,
    target.fdr = target_fdr,
    dist = cis_window,
    gene1.type = request$dataset$exposure_type,
    SNP.ref = snp_ref,
    n.cores = 1L,
    verbose = FALSE
  )
  artifact <- file.path(job_dir, sprintf("dandelion-chunk-%04d.rds", chunk_index))
  saveRDS(result, artifact, compress = "xz")
  artifact_paths <- c(artifact_paths, basename(artifact))

  for (exposure in result$gene1) {
    p_values <- result$mat.p[, exposure]
    valid <- !is.na(p_values)
    if (!any(valid)) next
    q_values <- rep(NA_real_, length(p_values))
    q_values[valid] <- stats::p.adjust(p_values[valid], method = "BH")
    package_significant <- result$mat.sig[, exposure] != 0
    derived_significant <- !is.na(q_values) & q_values <= target_fdr
    if (!identical(unname(as.logical(package_significant)), unname(as.logical(derived_significant)))) {
      stop(sprintf("BH parity validation failed for exposure %s.", exposure))
    }
    tested_pairs <- tested_pairs + sum(valid)
    keep <- which(derived_significant)
    nonsig <- which(valid & !derived_significant)
    if (length(nonsig) > 0L) keep <- unique(c(keep, nonsig[order(q_values[nonsig], p_values[nonsig])][seq_len(min(20L, length(nonsig)))]))
    for (row_index in keep) {
      candidate <- rownames(result$mat.p)[row_index]
      source_gene <- exposure
      source_type <- if (request$dataset$exposure_type == "SNP") "variant" else "gene"
      if (request$dataset$exposure_type == "SNP" && !is.null(snp_map)) {
        mapped <- snp_map$GeneSymbol[match(exposure, snp_map$SNP)]
        if (length(mapped) && !is.na(mapped) && nzchar(mapped)) {
          source_gene <- mapped
          source_type <- "gene_mapped_from_variant"
        }
      }
      record_index <- record_index + 1L
      all_records[[record_index]] <- list(
        exposure = exposure,
        exposure_type = request$dataset$exposure_type,
        source_gene = source_gene,
        source_node_type = source_type,
        candidate_gene = candidate,
        trans_p_value = as.numeric(p_trans[candidate, exposure]),
        gene_association_p_value = as.numeric(p_wes[candidate]),
        p_value = as.numeric(p_values[row_index]),
        q_value = as.numeric(q_values[row_index]),
        significant = isTRUE(derived_significant[row_index]),
        gene_association_threshold_pass = isTRUE(p_wes[candidate] <= gene_threshold)
      )
    }
  }
}

write_progress(96, "normalizing_results")
output <- list(
  contract_version = "1.0",
  method = "dandelion",
  package = list(name = "DANDELION", version = "0.1.0", license = "GPL-3.0"),
  multiple_testing = "Benjamini-Hochberg within each exposure",
  tested_pair_count = tested_pairs,
  returned_record_count = length(all_records),
  records = all_records,
  artifacts = artifact_paths,
  limitations = c(
    "Non-significant normalized output is limited to the top 20 pairs per exposure; complete package objects are retained as RDS artifacts.",
    "Candidate prioritisation is hypothesis-generating and does not establish causality."
  ),
  warnings = if (trans_spec$format == "rds") c("RDS input was loaded as one object and could not be streamed by exposure.") else character()
)
jsonlite::write_json(output, result_path, auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = TRUE)
write_progress(100, "completed")
