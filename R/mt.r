####################################################################################
# # # # Cannot do broad import of multtest. Its plot and summary export conflicts
# # # # conflicts with phylobase
# # # # @import multtest
####################################################################################
#' Wrapper of mt.maxT and mt.minP for the phyloseq package.
#'
#' @usage mt(X, classlabel, minPmaxT="minP", ...)
#'
#' @param X (Required). An \code{otuTable} or an object that contains an \code{otuTable}.
#'  In this multiple testing framework, different taxa correspond to variables
#'  (hypotheses), and samples to observations.
#'
#' @param classlabel (Required). Either a single character string indicating the column
#'  in the \code{sampleMap} of X, or, alternatively, an integer vector
#'  of length equal to the number of samples of X. The latter option is
#'  equivalent to \code{classlabel} of \code{\link[multtest]{mt.maxT}}. A 
#'  third option is for classlabel to be a 2-level factor, with length
#'  equal to the number of samples of X. If \code{classlabel} is a non-integer
#'  numeric class, \code{mt()} will attempt to coerce it to integer.
#'
#' @param minPmaxT (Optional). Whether to use \code{mt.minP} or \code{mt.maxT}.
#'  Default is to use \code{mt.minP}
#' 
#' @param ... (Optional). Additional arguments, forwarded to \code{mt.maxT} 
#' 
#' @return A dataframe with components specified in the documentation for
#'  \code{mt.maxT} or \code{mt.minP} respectively.
#'
#' @seealso \code{\link[multtest]{mt.maxT}}, \code{\link[multtest]{mt.minP}}
#'
#' @rdname mt-methods
#' @docType methods
#' @export
#'
#' @importFrom multtest mt.maxT
#' @importFrom multtest mt.minP
#'
#' @examples #
#' ## data(ex1)
#' ## mt(ex1, "Diet")
#' ## # Alternatively
#' ## dietfac <- factor(data.frame(sampleMap(ex1))[, "Diet"])
#' ## mt(ex1, dietfac)
setGeneric("mt", function(X, classlabel, minPmaxT="minP", ...) standardGeneric("mt") )
################################################################################
#' @aliases mt,otuTable,integer-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "integer"), function(X, classlabel, minPmaxT="minP", ...){
	# Guarantee proper orientation of abundance table, and coerce to matrix.
	if( !speciesAreRows(X) ){ X <- t(X)	}
	mt.phyloseq.internal(as(X, "matrix"), classlabel, minPmaxT, ...)
})
################################################################################
# Force numeric classlabel to be integer, pass-on
#' @aliases mt,otuTable,numeric-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "numeric"), function(X, classlabel, minPmaxT="minP", ...){
	mt(X, as(classlabel, "integer"), minPmaxT="minP", ...)
})
################################################################################
#' @aliases mt,phyloseq,ANY-method
#' @rdname mt-methods
setMethod("mt", c("phyloseq", "ANY"), function(X, classlabel, minPmaxT="minP", ...){
	if( !is.null(access(X, "sampleMap")) & class(classlabel)=="character" ){
		rawFactor  <- as(sampleMap(X), "data.frame")[, classlabel[1]]
		if( class(rawFactor) != "factor" ){
			rawFactor <- factor(rawFactor)
		}
		classlabel <- rawFactor
	}
	mt(otuTable(X), classlabel, minPmaxT, ...)
})
################################################################################
#' @aliases mt,otuTable,factor-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "factor"), function(X, classlabel, minPmaxT="minP", ...){
	if( length(levels(classlabel)) != 2 ){
		stop("classlabel argument should be (or specify in sampleMap) a two-level factor.\n",
		"In this case, the classlabel had ", length(levels(classlabel)), " levels (after coercion)."
		)
	}
	# integerize classlabel to logical suitable for mtmaxT
	classlabel <- (0:1)[classlabel]
	# Use mt dispatch with classlabel now a suitable classlabel
	mt(X, classlabel, minPmaxT, ...)
})
####################################################################################
# Internal function
# @aliases mt,matrix,integer-method
# not exported
# @keywords internal
mt.phyloseq.internal <- function(X, classlabel, minPmaxT="minP", ...){
	# require(multtest)
	if( minPmaxT == "minP" ){
		mt.minP(X, classlabel, ...)
	} else if( minPmaxT == "maxT" ){
		mt.maxT(X, classlabel, ...)	
	} else {
		print("Nothing calculated. minPmaxT argument must be either minP or maxT.")
	}
}
####################################################################################