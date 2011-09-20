####################################################################################
#' Wrapper of mt.maxT and mt.minP for the phyloseq package.
#'
#' @param X An \code{otuTable} or an object that contains an \code{otuTable}.
#'  In this multiple testing framework, different taxa correspond to variables
#'  (hypotheses), and samples to observations.
#'
#' @param classlabel Either a single character string indicating the column
#'  in the \code{sampleMap} of X, or, alternatively, an integer vector
#'  of length equal to the number of samples of X. The latter option is
#'  equivalent to \code{classlabel} of \code{mt.maxT} in \code{multtest}. A 
#'  third option is for classlabel to be a 2-level factor, with length
#'  equal to the number of samples of X. 
#'
#' @param ... Additional arguments, forwarded to \code{mt.maxT} 
#' 
#' @return A dataframe with components specified in the documentation for
#'  \code{mt.maxT} or \code{mt.minP} respectively.
#'
#' @seealso mt.maxT mt.minP
#'
#' @rdname mt-methods
#' @docType methods
#' @export
#'
#' @import multtest
#' @examples #
#' ## data(ex1)
#' ## mt(ex1, "Diet")
#' ## mt(otuSamTax(ex1), "Diet", minPmaxT="minsldzxc")
setGeneric("mt", function(X, classlabel, minPmaxT="minP", ...)
	standardGeneric("mt")
)
####################################################################################
# otuTree doesn't have a sampleMap. just grab otuTable hope 4 best.
#' @aliases mt,matrix,integer-method
#' @rdname mt-methods
setMethod("mt", c("matrix", "integer"), function(X, classlabel, minPmaxT="minP", ...){
	# require(multtest)
	if( minPmaxT == "minP" ){
		mt.minP(X, classlabel, minPmaxT, ...)
	} else if( minPmaxT == "maxT" ){
		mt.maxT(X, classlabel, minPmaxT, ...)	
	} else {
		print("Nothing calculated. minPmaxT argument must be either minP or maxT.")
	}
})
################################################################################
#' @aliases mt,ANY,integer-method
#' @rdname mt-methods
setMethod("mt", c("ANY", "integer"), function(X, classlabel, minPmaxT="minP", ...){
	X <- otuTable(X)
	if( !speciesAreRows(X) ){ X <- t(X)	}
	mt(as(X, "matrix"), classlabel, minPmaxT, ...)
})
################################################################################
#' @aliases mt,otuSam,character-method
#' @rdname mt-methods
setMethod("mt", c("otuSam", "character"), function(X, classlabel, minPmaxT="minP", ...){
	rawFactor <- as(sampleMap(X), "data.frame")[, classlabel]
	mt(otuTable(X), factor(rawFactor), minPmaxT, ...)
})
################################################################################
#' @aliases mt,ANY,factor-method
#' @rdname mt-methods
setMethod("mt", c("ANY", "factor"), function(X, classlabel, minPmaxT="minP", ...){
	# integerize classlabel to logical suitable for mtmaxT
	classlabel <- (0:1)[classlabel]
	# Use dispatch of mtmaxT with classlabel now a suitable classlabel
	mt(otuTable(X), classlabel, minPmaxT, ...)
})
################################################################################
