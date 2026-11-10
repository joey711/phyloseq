.asv_ids <- function(ids, argument="sequence IDs") {
  if (!is.character(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Expected unique, non-missing, non-empty ", argument, ".")
  }
}

.asv_dna <- function(sequences, unambiguous=TRUE) {
  if (!inherits(sequences, "XStringSet") && !is.character(sequences)) {
    stop("Expected named nucleotide sequences or an XStringSet.")
  }
  ids <- names(sequences)
  .asv_ids(ids)
  text <- toupper(as.character(sequences))
  text <- chartr("U", "T", text)
  alphabet <- if (unambiguous) "^[ACGT]+$" else "^[ACGTRYSWKMBDHVN]+$"
  if (anyNA(text) || any(!grepl(alphabet, text))) {
    stop(if (unambiguous) "ASVs must contain non-empty, unambiguous A/C/G/T sequences (U is converted to T)." else
      "Expected non-empty, unaligned IUPAC nucleotide sequences (U is converted to T).")
  }
  Biostrings::DNAStringSet(stats::setNames(text, ids))
}

.asv_scalar <- function(x, argument, minimum=0, maximum=Inf, whole=FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x > maximum || (whole && x != floor(x))) {
    stop("Invalid `", argument, "`: expected one finite ", if (whole) "whole " else "",
         "number in [", minimum, ", ", maximum, "].")
  }
}

#' Extract observed ASV sequences with strict feature-ID matching
#'
#' Extracts phyloseq refseq or TreeSummarizedExperiment referenceSeq. Sequences
#' must cover exactly the features in the input; no feature is silently dropped.
#' Sequence case is normalized and RNA U becomes DNA T. Ambiguous or gapped
#' observed sequences are rejected. Repeated sequences with distinct feature IDs
#' are retained and never merged implicitly.
#'
#' @param x A phyloseq, TreeSummarizedExperiment, named character vector or
#'   XStringSet containing the observed nucleotide sequences.
#' @param ids Optional character feature IDs in the desired output order.
#' @return A named DNAStringSet. Original input objects are unchanged.
#' @export
#' @examples
#' sequences <- c("01"="acgu", asv2="TTAA")
#' asv_sequences(sequences)
asv_sequences <- function(x, ids=NULL) {
  expected <- NULL
  if (inherits(x, "phyloseq")) {
    expected <- taxa_names(x)
    sequences <- refseq(x, errorIfNULL=FALSE)
  } else if (inherits(x, "TreeSummarizedExperiment")) {
    .composition_require("TreeSummarizedExperiment")
    expected <- rownames(x)
    sequences <- TreeSummarizedExperiment::referenceSeq(x)
  } else {
    sequences <- x
  }
  if (is.null(sequences)) stop("Observed ASV sequences are missing; supply refseq/referenceSeq first.")
  sequences <- .asv_dna(sequences)
  if (!is.null(expected)) {
    .asv_ids(expected, "feature IDs")
    if (!setequal(expected, names(sequences))) stop("Sequence IDs do not match all input feature IDs.")
    sequences <- sequences[expected]
  }
  if (!is.null(ids)) {
    .asv_ids(ids, "requested feature IDs")
    if (!all(ids %in% names(sequences))) stop("Requested feature IDs are absent from the input.")
    sequences <- sequences[ids]
  }
  sequences
}

#' Select ASVs by feature ID or observed nucleotide sequence
#'
#' Specify exactly one selector. Sequence selection retains every matching
#' feature, including different IDs carrying the same DNA. It never substitutes
#' a reference gene for an observed ASV. Container selection retains input feature
#' order, samples, annotations and existing transformation records, so subsequent
#' check_composition calls can detect assays made stale by feature removal.
#'
#' @inheritParams asv_sequences
#' @param sequences Character DNA/RNA sequences to match exactly after case and
#'   U-to-T normalization. All requested sequences must match at least one ASV.
#' @param both_strands If TRUE, sequence selection also matches reverse complements.
#' @param prune If FALSE, return a logical selection named by input feature ID.
#' @return With prune=TRUE, the selected input container, or a DNAStringSet for
#'   sequence-vector inputs. An empty selection is an error. phyloseq tree pruning
#'   follows the native prune_taxa singleton-tree behavior.
#' @export
#' @examples
#' sequences <- c("01"="ACGT", asv2="TTAA", asv3="ACGT")
#' select_asvs(sequences, sequences="ACGT")
select_asvs <- function(x, ids=NULL, sequences=NULL, both_strands=FALSE, prune=TRUE) {
  if (is.null(ids) == is.null(sequences)) stop("Supply exactly one of `ids` and `sequences`.")
  for (argument in c("both_strands", "prune")) {
    value <- get(argument)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("`", argument, "` must be TRUE or FALSE.")
  }
  observed <- asv_sequences(x)
  if (!is.null(ids)) {
    .asv_ids(ids, "requested feature IDs")
    if (!all(ids %in% names(observed))) stop("Requested feature IDs are absent from the input.")
    keep <- names(observed) %in% ids
  } else {
    if (!is.character(sequences) || !length(sequences)) stop("`sequences` must contain nucleotide strings.")
    names(sequences) <- paste0("selector", seq_along(sequences))
    requested <- .asv_dna(sequences)
    requested_text <- as.character(requested)
    observed_text <- as.character(observed)
    matched <- requested_text %in% observed_text
    if (both_strands) matched <- matched | as.character(Biostrings::reverseComplement(requested)) %in% observed_text
    if (!all(matched)) stop("Some requested sequences do not match any observed ASV.")
    if (both_strands) requested_text <- c(requested_text, as.character(Biostrings::reverseComplement(requested)))
    keep <- observed_text %in% requested_text
  }
  names(keep) <- names(observed)
  if (!prune) return(keep)
  if (!any(keep)) stop("No ASVs selected.")
  if (inherits(x, "phyloseq")) return(prune_taxa(keep, x))
  if (inherits(x, "TreeSummarizedExperiment")) return(x[keep, , drop=FALSE])
  observed[keep]
}
