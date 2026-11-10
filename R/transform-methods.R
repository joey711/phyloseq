#' Thresholded rank transformation.
#'
#' The lowest \code{thresh} values in \code{x} all get the value 'thresh'.
#'
#' @usage threshrank(x, thresh, keep0s=FALSE, ...)
#'
#' @param x (Required). Numeric vector to transform.
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#'
#' @return A ranked, (optionally) thresholded numeric vector with length equal to
#'  \code{x}. Default arguments to \code{rank} are used, unless provided as
#'  additional arguments.
#'
#' @seealso \code{\link{transform_sample_counts}}, \code{\link{rank}}, \code{\link{threshrankfun}}
#' @export
#' @examples #
#' (a_vector <- sample(0:10, 100, TRUE))
#' threshrank(a_vector, 5, keep0s=TRUE)
#' data(GlobalPatterns)
#' GP <- GlobalPatterns
#' ## These three approaches result in identical otu_table
#' (x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
#' (x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
#' identical(x1, x2)
#' (x3 <- otu_table(apply(otu_table(GP), 2, threshrank, thresh=500), taxa_are_rows(GP)) )
#' identical(x1, x3)
threshrank <- function(x, thresh, keep0s=FALSE, ...){
	if( keep0s ){ index0 <- which(x == 0) }
	x <- rank(x, ...)
	thresh <- thresh[1]
	x[x<thresh] <- thresh
	if( keep0s ){ x[index0] <- 0 }
	return(x)
}
####################################################################################
#' A closure version of the \code{threshrank} function.
#'
#' Takes the same arguments as \code{\link{threshrank}}, except for \code{x},
#' because the output is a single-argument function rather than a rank-transformed numeric.
#' This is useful for higher-order functions that require a single-argument function as input,
#' like \code{\link{transform_sample_counts}}.
#'
#' @usage threshrankfun(thresh, keep0s=FALSE, ...)
#'
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#'
#' @return A single-argument function with the options to \code{\link{threshrank}} set.
#'
#' @seealso \code{\link{transform_sample_counts}}, \code{\link{threshrankfun}},
#'  \code{\link{threshrank}}
#' @export
#' @examples
#' data(esophagus)
#' x1 = transform_sample_counts(esophagus, threshrankfun(50))
#' otu_table(x1)
#' x2 = transform_sample_counts(esophagus, rank)
#' otu_table(x2)
#' identical(x1, x2)
threshrankfun <- function(thresh, keep0s=FALSE, ...){
	function(x){
		threshrank(x, thresh, keep0s=keep0s, ...)
	}
}
################################################################################
#' Transpose \code{\link{otu_table-class}} or \code{\link{phyloseq-class}}
#'
#' Extends the base transpose method, \code{\link[base]{t}}.
#'
#' @usage t(x)
#'
#' @param x An \code{otu_table} or \code{\link{phyloseq-class}}.
#'
#' @return The class of the object returned by \code{t} matches
#' the class of the argument, \code{x}. The \code{otu_table} is
#' transposed, and \code{\link{taxa_are_rows}} value is toggled.
#'
#' @name t
#' @rdname transpose-methods
#' @docType methods
#' @export
#' @examples
#' data(GlobalPatterns)
#' otu_table(GlobalPatterns)
#' t( otu_table(GlobalPatterns) )
setGeneric("t")
#' @aliases t,otu_table-method
#' @rdname transpose-methods
setMethod("t", signature("otu_table"), function(x){
	#new("otu_table", t(x@.Data), taxa_are_rows = (!taxa_are_rows(x)))
	x <- otu_table( t(as(x, "matrix")), taxa_are_rows=(!taxa_are_rows(x)) )
	return(x)
})
################################################################################
#' @aliases t,phyloseq-method
#' @rdname transpose-methods
setMethod("t", signature("phyloseq"), function(x){
	x@otu_table <- t( otu_table(x) )
	return(x)
})
################################################################################
#' Transform abundance data in an \code{otu_table}, sample-by-sample.
#'
#' This function transforms the sample counts of a taxa
#' abundance matrix according to a user-provided function.
#' The counts of each sample will be transformed individually. No sample-sample
#' interaction/comparison is possible by this method.
#'
#' @usage transform_sample_counts(physeq, fun, ...)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} of \code{\link{otu_table-class}}.
#'
#' @param fun (Required). A single-argument function that will be applied
#'  to the abundance counts of each sample.
#'  Can be an anonymous \code{\link[base]{function}}.
#'  It must return a numeric vector of the same length as its input.
#'  The function is called once per sample, using the actual sample counts
#'  and any additional arguments in \code{...}.
#'
#' @param ... (Optional). Additional, optionally-named, arguments passed to
#'  \code{fun} during transformation of abundance data.
#'
#' @return A transformed \code{otu_table} -- or \code{phyloseq} object with its
#'  transformed \code{otu_table}.
#'  Orientation, dimensions, and sample/taxon names are preserved, as are
#'  the other components of a \code{phyloseq} object. Results are assigned
#'  by position; names returned by \code{fun} do not change taxon identities.
#'  Filtering/trimming can follow, for which the
#'  \code{\link{genefilter_sample}} and \code{\link{prune_taxa}} functions
#'  are suggested.
#'
#' @seealso \code{\link{threshrankfun}}, \code{\link{rank}}, \code{\link{log}}
#'
#' @docType methods
#' @aliases transform_sample_counts transformSampleCounts
#' @rdname transformcounts
#' @export
#'
#' @examples #
#' data(esophagus)
#' x1 = transform_sample_counts(esophagus, threshrankfun(50))
#' head(otu_table(x1), 10)
#' x2 = transform_sample_counts(esophagus, rank)
#' head(otu_table(x2), 10)
#' identical(x1, x2)
#' x3 = otu_table(esophagus) + 5
#' x3 = transform_sample_counts(x3, log)
#' head(otu_table(x3), 10)
#' x4 = transform_sample_counts(esophagus, function(x) round(x^2.2, 0))
#' head(otu_table(x4), 10)
transform_sample_counts <- function(physeq, fun, ...){
  if (!is.function(fun)) {
    stop("`fun` must be a function.")
  }

  # Work in one orientation, then restore the original layout on output.
  taxa_rows <- taxa_are_rows(physeq)
  counts <- as(otu_table(physeq), "matrix")
  if (!taxa_rows) counts <- t(counts)

  transformed <- lapply(seq_len(ncol(counts)), function(i) {
    values <- fun(counts[, i], ...)
    if (!is.numeric(values) || !is.null(dim(values)) || length(values) != nrow(counts)) {
      stop("`fun` must return a numeric vector of the same length as its input.")
    }
    values
  })

  # Explicit dimensions also cover a single taxon or a single sample.
  counts <- matrix(unlist(transformed, use.names=FALSE),
                   nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts))
  if (!taxa_rows) counts <- t(counts)
  otu_table(physeq) <- otu_table(counts, taxa_are_rows=taxa_rows)
  physeq
}
####################################################################################
#' @rdname transformcounts
#' @export
transformSampleCounts <- transform_sample_counts
####################################################################################
############################################################
