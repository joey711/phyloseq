# Internal extraction uses taxa in rows and samples in columns throughout.
.composition_matrix <- function(x, assay="counts") {
  .composition_name(assay, "assay")
  if (inherits(x, "phyloseq") || inherits(x, "otu_table")) {
    values <- as(otu_table(x), "matrix")
    if (!taxa_are_rows(x)) values <- t(values)
  } else if (inherits(x, "SummarizedExperiment")) {
    .composition_require("SummarizedExperiment")
    if (!assay %in% SummarizedExperiment::assayNames(x)) {
      stop("Assay not found: ", assay)
    }
    values <- as.matrix(SummarizedExperiment::assay(x, assay))
  } else if (is.matrix(x)) {
    values <- x
  } else {
    stop("Expected phyloseq, otu_table, SummarizedExperiment, or a taxa-by-sample matrix.")
  }
  if (!is.numeric(values) || is.complex(values) || any(dim(values) == 0L)) {
    stop("Composition data must be a non-empty, real numeric matrix.")
  }
  values
}

.composition_require <- function(package) {
  if (!requireNamespace(package, quietly=TRUE)) {
    stop("This operation requires the optional package '", package,
         "'. Install it with BiocManager::install(\"", package, "\").")
  }
}

.composition_name <- function(value, argument) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("`", argument, "` must be one non-empty character string.")
  }
}

.composition_id_issues <- function(values) {
  issues <- character()
  for (margin in c("taxon", "sample")) {
    ids <- if (margin == "taxon") rownames(values) else colnames(values)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      issues <- c(issues, paste("Expected unique, non-missing, non-empty", margin, "IDs."))
    }
  }
  issues
}

.composition_fingerprints <- function(values, version=2L) {
  hashes <- vapply(seq_len(ncol(values)), function(i) {
    column <- values[, i]
    # Taxon IDs are checked separately and ordered before verification.
    # Version 1 retains names to verify assays created by the initial release.
    if (version == 2L) column <- unname(column)
    digest::digest(column, algo="xxhash64")
  }, character(1))
  stats::setNames(hashes, colnames(values))
}

.composition_provenance_issues <- function(x, assay, values, seen=character()) {
  if (!inherits(x, "SummarizedExperiment")) return(character())
  info <- S4Vectors::metadata(x)$phyloseqV2
  if (!is.list(info) || !is.list(info$transformations)) return(character())
  record <- info$transformations[[assay]]
  if (is.null(record)) return(character())
  if (assay %in% seen) return("Transformation provenance contains a cycle; recompute this assay.")
  if (!is.list(record) || !is.character(record$source_assay) ||
      length(record$source_assay) != 1L || is.na(record$source_assay) ||
      is.null(record$taxa) || is.null(record$source_hashes) || is.null(record$result_hashes)) {
    return("Transformation provenance is incomplete; recompute this assay.")
  }
  if (!setequal(rownames(values), record$taxa)) {
    return("The taxon set changed after transformation; recompute this assay from the source data.")
  }
  if (!record$source_assay %in% SummarizedExperiment::assayNames(x)) {
    return("The source assay is missing; transformation provenance cannot be verified.")
  }
  version <- record$fingerprint_version
  if (is.null(version)) version <- 1L
  if (!is.numeric(version) || length(version) != 1L || is.na(version) ||
      !version %in% c(1L, 2L)) {
    return("Unsupported transformation fingerprint version; recompute this assay.")
  }
  source <- .composition_matrix(x, record$source_assay)
  result <- values
  if (!identical(rownames(values), record$taxa)) {
    source <- source[record$taxa, , drop=FALSE]
    result <- values[record$taxa, , drop=FALSE]
  }
  expected_source <- record$source_hashes[colnames(values)]
  expected_result <- record$result_hashes[colnames(values)]
  if (anyNA(expected_source) || anyNA(expected_result) ||
      !identical(unname(.composition_fingerprints(source, version)), unname(expected_source)) ||
      !identical(unname(.composition_fingerprints(result, version)), unname(expected_result))) {
    return("Source or transformed values changed after transformation; recompute this assay.")
  }
  .composition_provenance_issues(x, record$source_assay, source, seen=c(seen, assay))
}

#' Check abundance, proportion or CLR data before compositional analysis
#'
#' Checks values and identifiers without silently removing taxa or samples.
#' These checks validate numerical requirements, not the biological suitability
#' of a normalization or zero-replacement strategy. A zero-sum vector alone does
#' not prove that data were obtained by CLR transformation.
#' For assays created by transform_composition, recorded fingerprints also
#' detect changed source/output values or a changed taxon set. Sample subsets
#' and permutations are supported. External assays without such records receive
#' numerical checks only. The Observed and Zeros fields count positive and zero
#' entries; on the CLR scale these do not describe taxon presence or absence.
#'
#' @param x A phyloseq, otu_table, SummarizedExperiment (including
#'   TreeSummarizedExperiment), or numeric matrix with taxa in rows and samples
#'   in columns. Matrices must have unique row and column names.
#' @param assay Name of the assay to read for SummarizedExperiment inputs.
#' @param scale Expected input scale. Abundances must be finite, nonnegative and
#'   have positive mass in every sample. Counts must additionally be whole
#'   numbers. Proportions must sum to one per sample. CLR coordinates may be
#'   negative, must be centered per sample, and require at least two taxa.
#' @param tolerance Finite positive absolute tolerance for whole-number,
#'   proportion-sum and CLR-mean checks.
#' @param stop_on_error If TRUE, stop with the validation issues.
#' @return A list with valid, scale, issues, per-sample diagnostics in samples,
#'   and IDs of all-zero taxa in zero_taxa. All-zero taxa are reported, not removed.
#' @export
#' @examples
#' data(esophagus)
#' report <- check_composition(esophagus, scale="counts")
#' report$valid
#' head(report$samples)
check_composition <- function(x, assay="counts",
                              scale=c("abundance", "counts", "proportion", "clr"),
                              tolerance=1e-8, stop_on_error=FALSE) {
  scale <- match.arg(scale)
  if (!is.numeric(tolerance) || is.complex(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) stop("`tolerance` must be finite and positive.")
  if (!is.logical(stop_on_error) || length(stop_on_error) != 1L || is.na(stop_on_error)) {
    stop("`stop_on_error` must be TRUE or FALSE.")
  }
  values <- .composition_matrix(x, assay)
  issues <- c(.composition_id_issues(values), .composition_provenance_issues(x, assay, values))
  finite <- is.finite(values)
  nonfinite <- colSums(!finite)
  safe <- values
  if (any(nonfinite > 0)) safe[!finite] <- 0
  totals <- colSums(safe)
  centers <- colMeans(safe)
  totals[nonfinite > 0] <- centers[nonfinite > 0] <- NA_real_
  negative <- colSums(safe < 0)
  is_zero <- finite & values == 0
  zeros <- colSums(is_zero)
  observed <- nrow(values) - nonfinite - negative - zeros
  noninteger <- if (is.integer(safe)) {
    stats::setNames(rep(0, ncol(safe)), colnames(safe))
  } else {
    colSums(abs(safe - round(safe)) > tolerance)
  }
  valid <- nonfinite == 0
  if (any(!valid)) issues <- c(issues, "Non-finite values (NA, NaN or Inf) are present.")
  if (scale == "clr") {
    centered <- !is.na(centers) & abs(centers) <= tolerance
    valid <- valid & centered & nrow(values) >= 2L
    if (any(!centered)) issues <- c(issues, "CLR coordinates must have mean zero in each sample.")
    if (nrow(values) < 2L) issues <- c(issues, "CLR requires at least two taxa.")
  } else {
    valid <- valid & negative == 0 & observed > 0
    if (any(negative > 0)) issues <- c(issues, "Negative abundances are present.")
    if (any(observed == 0)) issues <- c(issues, "Samples with no positive abundance are present.")
    if (scale == "counts") {
      valid <- valid & noninteger == 0
      if (any(noninteger > 0)) issues <- c(issues, "Counts must be whole numbers.")
    }
    if (scale == "proportion") {
      closed <- !is.na(totals) & abs(totals - 1) <= tolerance
      valid <- valid & closed
      if (any(!closed)) issues <- c(issues, "Proportions must sum to one in each sample.")
    }
  }
  ids <- colnames(values)
  if (is.null(ids)) ids <- rep(NA_character_, ncol(values))
  report <- list(valid=length(issues) == 0L, scale=scale, issues=issues,
                 samples=data.frame(Sample=ids, Total=totals, Observed=observed,
                                    Zeros=zeros,
                                    ZeroFraction=zeros / nrow(values),
                                    NonFinite=nonfinite, Negative=negative,
                                    NonInteger=noninteger, Center=centers, Valid=valid,
                                    row.names=NULL, stringsAsFactors=FALSE),
                 zero_taxa=rownames(values)[rowSums(is_zero) == ncol(values)])
  if (stop_on_error && !report$valid) stop(paste(report$issues, collapse="\n"))
  report
}

#' Convert phyloseq data to a TreeSummarizedExperiment
#'
#' Keeps original values, identifiers, sample metadata, taxonomy, tree and
#' reference sequences. TreeSummarizedExperiment is an optional dependency.
#' This conversion does not infer whether the input contains counts or a
#' previously transformed abundance table; select an appropriate assay name.
#'
#' @param physeq A phyloseq or otu_table object.
#' @param assay Name for the original abundance assay.
#' @return A TreeSummarizedExperiment with taxa in rows and samples in columns.
#' @export
#' @examples
#' if (requireNamespace("TreeSummarizedExperiment", quietly=TRUE)) {
#'   data(esophagus)
#'   tse <- as_tree_experiment(esophagus)
#' }
as_tree_experiment <- function(physeq, assay="counts") {
  if (!inherits(physeq, "phyloseq") && !inherits(physeq, "otu_table")) {
    stop("`physeq` must be a phyloseq or otu_table object.")
  }
  values <- .composition_matrix(physeq, assay)
  issues <- .composition_id_issues(values)
  if (length(issues)) stop(paste(issues, collapse="\n"))
  .composition_require("TreeSummarizedExperiment")
  .composition_require("S4Vectors")
  taxonomy <- access(physeq, "tax_table")
  samples <- access(physeq, "sam_data")
  rows <- S4Vectors::DataFrame(row.names=rownames(values))
  columns <- S4Vectors::DataFrame(row.names=colnames(values))
  if (!is.null(taxonomy)) {
    if (anyDuplicated(rownames(taxonomy)) || !setequal(rownames(taxonomy), rownames(values))) {
      stop("Taxonomy taxon IDs do not match the abundance table.")
    }
    rows <- S4Vectors::DataFrame(as(taxonomy, "matrix")[rownames(values), , drop=FALSE])
  }
  if (!is.null(samples)) {
    if (anyDuplicated(rownames(samples)) || !setequal(rownames(samples), colnames(values))) {
      stop("Metadata sample IDs do not match the abundance table.")
    }
    columns <- S4Vectors::DataFrame(as(samples, "data.frame")[colnames(values), , drop=FALSE])
  }
  sequences <- access(physeq, "refseq")
  if (!is.null(sequences)) {
    if (anyDuplicated(names(sequences)) || !setequal(names(sequences), rownames(values))) {
      stop("Reference sequence taxon IDs do not match the abundance table.")
    }
    sequences <- sequences[rownames(values)]
  }
  tree <- access(physeq, "phy_tree")
  if (!is.null(tree) && (anyDuplicated(tree$tip.label) || !setequal(tree$tip.label, rownames(values)))) {
    stop("Tree tip IDs do not match the abundance table.")
  }
  result <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays=stats::setNames(list(values), assay), rowData=rows, colData=columns,
    rowTree=tree, referenceSeq=sequences)
  if (!identical(rownames(result), rownames(values)) ||
      !identical(colnames(result), colnames(values))) {
    stop("Conversion changed taxon or sample identifiers; refusing to return altered data.")
  }
  result
}

.close_composition <- function(values) {
  # Scaling first avoids overflowing the sample total for large finite inputs.
  scaled <- sweep(values, 2, apply(values, 2, max), "/")
  sweep(scaled, 2, colSums(scaled), "/")
}

#' Add a compositional transformation while retaining the source assay
#'
#' phyloseq and otu_table inputs are converted to TreeSummarizedExperiment;
#' SummarizedExperiment inputs retain their class and other assays. The source
#' assay is never overwritten. No samples or taxa are removed automatically.
#'
#' CLR uses natural logarithms, log(x) minus the sample mean of log(x).
#' Zero replacement is explicit: pseudocount adds the chosen value to every
#' component on the input scale; multiplicative first closes each sample,
#' replaces each zero by delta, and scales positive proportions by
#' (1 - number_of_zeros * delta), preserving ratios between positive components.
#' delta must be positive and smaller than 1 / number_of_zeros in every affected
#' sample. Neither strategy can distinguish structural from sampling zeros.
#' A fixed additive pseudocount breaks scale invariance. Sensitivity to zero
#' replacement should be assessed for the intended biological analysis.
#'
#' @param x A phyloseq, otu_table, or SummarizedExperiment object.
#' @param method relative (proportions), hellinger (square root of proportions),
#'   or clr (centered log ratio). Hellinger is not a log-ratio transformation.
#' @param assay Source assay name; also used for conversion from phyloseq.
#' @param name Name of the new assay. Existing assays cannot be overwritten.
#' @param zero_method CLR zero policy: error, pseudocount or multiplicative.
#' @param pseudocount Explicit finite positive scalar for pseudocount replacement.
#' @param delta Explicit proportion assigned to each zero for multiplicative replacement.
#' @param input_scale abundance, counts or proportion; checked before transformation.
#' @return A TreeSummarizedExperiment (for phyloseq/otu_table inputs), or the
#'   same SummarizedExperiment subclass as x, with an additional named assay.
#'   metadata(result)$phyloseqV2$transformations records the source and parameters.
#'   Per-sample source and result fingerprints support later validation with
#'   check_composition. Recompute after changing the taxon set or source values.
#' @references
#' \url{https://search.r-project.org/CRAN/refmans/compositions/html/clr.html}
#' \url{https://doi.org/10.1023/A:1023866030544}
#' @export
#' @examples
#' if (requireNamespace("TreeSummarizedExperiment", quietly=TRUE)) {
#'   data(esophagus)
#'   tse <- transform_composition(esophagus, "relative")
#'   tse <- transform_composition(tse, "clr", zero_method="pseudocount", pseudocount=1)
#' }
transform_composition <- function(x, method=c("relative", "hellinger", "clr"),
                                  assay="counts", name=method,
                                  zero_method=c("error", "pseudocount", "multiplicative"),
                                  pseudocount=NULL, delta=NULL,
                                  input_scale=c("abundance", "counts", "proportion")) {
  method <- match.arg(method)
  zero_method <- match.arg(zero_method)
  input_scale <- match.arg(input_scale)
  .composition_name(name, "name")
  check_composition(x, assay, scale=input_scale, stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  if (method != "clr" && (zero_method != "error" || !is.null(pseudocount) || !is.null(delta))) {
    stop("Zero replacement parameters apply only to CLR.")
  }
  if (zero_method != "pseudocount" && !is.null(pseudocount)) {
    stop("`pseudocount` requires zero_method='pseudocount'.")
  }
  if (zero_method != "multiplicative" && !is.null(delta)) {
    stop("`delta` requires zero_method='multiplicative'.")
  }
  if (method == "clr") {
    if (nrow(values) < 2L) stop("CLR requires at least two taxa.")
    if (zero_method == "error" && any(values == 0)) {
      stop("CLR requires positive data. Choose an explicit zero replacement strategy.")
    }
    if (zero_method == "pseudocount") {
      if (!is.numeric(pseudocount) || is.complex(pseudocount) || length(pseudocount) != 1L ||
          !is.finite(pseudocount) || pseudocount <= 0) {
        stop("`pseudocount` must be one finite positive number on the input scale.")
      }
      values <- values + pseudocount
    }
    if (zero_method == "multiplicative") {
      if (!is.numeric(delta) || is.complex(delta) || length(delta) != 1L ||
          !is.finite(delta) || delta <= 0 || delta >= 1) {
        stop("`delta` must be one finite number strictly between zero and one.")
      }
      zeros <- values == 0
      mass <- colSums(zeros) * delta
      if (any(mass >= 1)) stop("`delta` leaves no positive mass for observed components.")
      values <- sweep(.close_composition(values), 2, 1 - mass, "*")
      values[zeros] <- delta
    }
    if (any(!is.finite(values)) || any(values <= 0)) {
      stop("Zero replacement produced non-finite or non-positive values; review the parameters and dynamic range.")
    }
    logs <- log(values)
    transformed <- sweep(logs, 2, colMeans(logs), "-")
  } else {
    transformed <- .close_composition(values)
    if (method == "hellinger") transformed <- sqrt(transformed)
  }
  if (any(!is.finite(transformed))) stop("Transformation produced non-finite values.")
  if (inherits(x, "phyloseq") || inherits(x, "otu_table")) x <- as_tree_experiment(x, assay)
  if (!inherits(x, "SummarizedExperiment")) {
    stop("`x` must be phyloseq, otu_table, or SummarizedExperiment; matrices can be checked but are not assay containers.")
  }
  if (name %in% SummarizedExperiment::assayNames(x)) {
    stop("Refusing to overwrite existing assay: ", name, ". Choose another `name`.")
  }
  SummarizedExperiment::assay(x, name) <- transformed
  info <- S4Vectors::metadata(x)
  if (!is.null(info$phyloseqV2) && !is.list(info$phyloseqV2)) {
    stop("metadata(x)$phyloseqV2 must be a list to record transformations.")
  }
  if (is.null(info$phyloseqV2)) info$phyloseqV2 <- list()
  if (!is.null(info$phyloseqV2$transformations) && !is.list(info$phyloseqV2$transformations)) {
    stop("metadata(x)$phyloseqV2$transformations must be a list.")
  }
  info$phyloseqV2$transformations[[name]] <- list(
    method=method, source_assay=assay, input_scale=input_scale,
    zero_method=if (method == "clr") zero_method else NULL,
    pseudocount=pseudocount, delta=delta,
    log_base=if (method == "clr") "e" else NULL,
    taxa=rownames(values), samples=colnames(values), fingerprint_version=2L,
    source_hashes=.composition_fingerprints(.composition_matrix(x, assay)),
    result_hashes=.composition_fingerprints(transformed))
  S4Vectors::metadata(x) <- info
  x
}

#' Calculate Aitchison distances from explicitly selected CLR coordinates
#'
#' No zero replacement or normalization is performed implicitly. Euclidean
#' distance between sample CLR vectors equals their Aitchison distance. CLR
#' centering is checked, but numerical checks cannot establish transformation
#' provenance. Choose the intended assay explicitly when its name is not clr.
#'
#' @inheritParams check_composition
#' @return A dist object, with sample IDs as labels.
#' @export
aitchison_distance <- function(x, assay="clr", tolerance=1e-8) {
  check_composition(x, assay, scale="clr", tolerance=tolerance, stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  if (ncol(values) < 2L) stop("Aitchison distances require at least two samples.")
  result <- stats::dist(t(values))
  if (any(!is.finite(result))) stop("Aitchison distances are non-finite; review the coordinate magnitudes.")
  result
}
