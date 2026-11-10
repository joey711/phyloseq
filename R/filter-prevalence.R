#' Report taxon prevalence across samples
#'
#' A taxon is present in a sample when its abundance is strictly greater than
#' detection. Prevalence is the number of such samples divided by the total
#' number of samples, including all-zero samples. No samples are removed and
#' missing values are rejected rather than changing the denominator.
#'
#' Use original counts or nonnegative abundances before zero replacement.
#' Detection is expressed on the selected assay's scale: a read-count threshold
#' and a proportion threshold have different meanings. Recorded CLR assays are
#' rejected because CLR coordinates do not represent presence. For external
#' assays without transformation records, the caller must establish the scale
#' and suitability of the data.
#'
#' @inheritParams check_composition
#' @param detection A finite nonnegative number. Presence requires abundance
#'   strictly greater than this value; the default 0 counts positive entries.
#' @return A data frame in input taxon order, with Taxon (ID), Present (number
#'   of samples exceeding detection), Total (number of samples) and Prevalence
#'   (Present divided by Total). All-zero taxa and samples are retained.
#' @seealso \code{\link{filter_prevalence}}, \code{\link{check_composition}}
#' @export
#' @examples
#' data(esophagus)
#' head(taxa_prevalence(esophagus))
taxa_prevalence <- function(x, detection=0, assay="counts") {
  if (!is.numeric(detection) || is.complex(detection) || length(detection) != 1L ||
      !is.finite(detection) || detection < 0) {
    stop("`detection` must be one finite nonnegative number.")
  }
  values <- .composition_matrix(x, assay)
  issues <- .composition_id_issues(values)
  if (length(issues)) stop(paste(issues, collapse="\n"))
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Prevalence requires finite, nonnegative abundances; missing values are not allowed.")
  }
  if (inherits(x, "SummarizedExperiment")) {
    info <- S4Vectors::metadata(x)$phyloseqV2
    if (is.list(info) && is.list(info$transformations)) {
      record <- info$transformations[[assay]]
      if (is.list(record) && identical(record$method, "clr")) {
        stop("CLR coordinates cannot define taxon presence; use a nonnegative abundance assay.")
      }
    }
    issues <- .composition_provenance_issues(x, assay, values)
    if (length(issues)) stop(paste(issues, collapse="\n"))
  }
  present <- rowSums(values > detection)
  data.frame(Taxon=rownames(values), Present=unname(present), Total=ncol(values),
             Prevalence=unname(present) / ncol(values),
             row.names=NULL, stringsAsFactors=FALSE)
}

#' Filter taxa by prevalence without removing samples
#'
#' Retains taxa present in at least the requested fraction or number of samples.
#' Specify exactly one of min_prevalence and min_samples. Retention includes the
#' threshold boundary, whereas presence requires abundance strictly greater than
#' detection. The default detection=0 treats any positive abundance as present.
#' Prevalence is computed across all samples, including all-zero samples.
#'
#' Filtering only subsets taxa: it does not renormalize abundances, replace zeros,
#' remove samples or recompute derived assays. An empty sample may remain or be
#' created; use check_composition before subsequent transformations.
#' Filter before compositional transformations when possible. For existing
#' derived assays, transformation records are retained unchanged so that
#' check_composition detects the changed taxon set. Recompute those assays from
#' the filtered source data using a new output assay name.
#'
#' phyloseq and otu_table inputs use prune_taxa, including its handling of trees
#' with fewer than two retained taxa. SummarizedExperiment inputs use native row
#' subsetting, retaining all assays and associated metadata.
#'
#' @inheritParams taxa_prevalence
#' @param min_prevalence Minimum fraction of samples, between 0 and 1 inclusive.
#'   Supply this or min_samples, not both. There is no default retention threshold.
#' @param min_samples Minimum whole number of samples, between 0 and the number
#'   of samples inclusive. An explicit alternative to min_prevalence.
#' @param prune If TRUE (default), return the filtered input. If FALSE, return
#'   a named logical vector indicating which taxa meet the threshold.
#' @return With prune=TRUE, an object of the input class, retaining the original
#'   sample order, taxon order and table orientation. If no taxa remain, an error
#'   is raised. With prune=FALSE, a logical vector named by taxon IDs, including
#'   the case where every entry is FALSE. The input is never modified.
#' @seealso \code{\link{taxa_prevalence}}, \code{\link{prune_taxa}},
#'   \code{\link{transform_composition}}
#' @export
#' @examples
#' data(esophagus)
#' # Illustrative thresholds; choose them for the study and input scale.
#' filtered <- filter_prevalence(esophagus, min_prevalence=0.5)
#' keep <- filter_prevalence(esophagus, min_samples=2, prune=FALSE)
#' table(keep)
filter_prevalence <- function(x, min_prevalence=NULL, min_samples=NULL,
                              detection=0, assay="counts", prune=TRUE) {
  if (is.null(min_prevalence) == is.null(min_samples)) {
    stop("Supply exactly one of `min_prevalence` and `min_samples`.")
  }
  if (!is.logical(prune) || length(prune) != 1L || is.na(prune)) {
    stop("`prune` must be TRUE or FALSE.")
  }
  if (!is.null(min_prevalence) &&
      (!is.numeric(min_prevalence) || is.complex(min_prevalence) ||
       length(min_prevalence) != 1L || !is.finite(min_prevalence) ||
       min_prevalence < 0 || min_prevalence > 1)) {
    stop("`min_prevalence` must be one finite number between 0 and 1.")
  }
  if (!is.null(min_samples) &&
      (!is.numeric(min_samples) || is.complex(min_samples) ||
       length(min_samples) != 1L || !is.finite(min_samples) ||
       min_samples < 0 || min_samples != floor(min_samples))) {
    stop("`min_samples` must be one finite nonnegative whole number.")
  }
  report <- taxa_prevalence(x, detection=detection, assay=assay)
  if (!is.null(min_samples) && min_samples > report$Total[1L]) {
    stop("`min_samples` cannot exceed the number of samples.")
  }
  # Compare proportions directly: ceiling(p * n) can round up at exact boundaries.
  keep <- if (is.null(min_samples)) {
    report$Prevalence >= min_prevalence
  } else {
    report$Present >= min_samples
  }
  names(keep) <- report$Taxon
  if (!prune) return(keep)
  if (!any(keep)) {
    stop("No taxa meet the prevalence threshold; use `prune=FALSE` to inspect the selection.")
  }
  if (all(keep)) return(x)
  if (inherits(x, "phyloseq") || inherits(x, "otu_table")) {
    return(prune_taxa(keep, x))
  }
  x[keep, , drop=FALSE]
}
