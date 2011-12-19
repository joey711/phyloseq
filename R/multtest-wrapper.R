####################################################################################
# # # # Avoiding full import of multtest to mitigate potential conflicts
####################################################################################
#' Multiple testing of taxa abundancesa acccording to a sample variate
#'
#' @usage mt(physeq, classlabel, minPmaxT="minP", ...)
#'
#' @param physeq (Required). \code{\link{otuTable-class}} or \code{\link{phyloseq-class}}.
#'  In this multiple testing framework, different taxa correspond to variables
#'  (hypotheses), and samples to observations.
#'
#' @param classlabel (Required). A single character index of the sample-variable
#'  in the \code{\link{sampleData}} of \code{physeq} that will be used for multiple testing. 
#'  Alternatively, \code{classlabel} can be a custom integer, character, or factor with
#'  length equal to \code{nsamples(physeq)}. In either scenario -- a sample variable
#'  within the phyloseq-class object or a custom vector -- \code{classlabel} must
#'  be coercable to a logical (i.e. two-level factor). A variable with more than
#'  two states is not appropriate for this category of test.
#'
#' @param minPmaxT (Optional). Character string. \code{"mt.minP"} or \code{"mt.maxT"}.
#'  Default is to use \code{\link[multtest]{mt.minP}}.
#' 
#' @param ... (Optional). Additional arguments, forwarded to
#'  \code{\link[multtest]{mt.maxT}} or \code{\link[multtest]{mt.minP}}
#' 
#' @return A dataframe with components specified in the documentation for
#'  \code{\link[multtest]{mt.maxT}} or \code{\link[multtest]{mt.minP}}, respectively.
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
#' ## dietfac <- factor(data.frame(sampleData(ex1))[, "Diet"])
#' ## mt(ex1, dietfac)
setGeneric("mt", function(physeq, classlabel, minPmaxT="minP", ...) standardGeneric("mt") )
################################################################################
#' @aliases mt,phyloseq,ANY-method
#' @rdname mt-methods
setMethod("mt", c("phyloseq", "ANY"), function(physeq, classlabel, minPmaxT="minP", ...){
	# If sampleData slot is non-empty, and the classlabel is a character-class
	# length(classlabel) == 1
	if( !is.null(sampleData(physeq, FALSE)) & class(classlabel)=="character" & length(classlabel)==1 ){
		rawFactor  <- as(sampleData(physeq), "data.frame")[, classlabel[1]]
		if( class(rawFactor) != "factor" ){
			rawFactor <- factor(rawFactor)
		}
		classlabel <- rawFactor
	} # Either way, dispatch on otuTable(physeq)
	mt(otuTable(physeq), classlabel, minPmaxT, ...)
})
################################################################################
#' @aliases mt,otuTable,integer-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "integer"), function(physeq, classlabel, minPmaxT="minP", ...){
	# Guarantee proper orientation of abundance table, and coerce to matrix.
	if( !speciesAreRows(physeq) ){ physeq <- t(physeq)	}
	mt.phyloseq.internal(as(physeq, "matrix"), classlabel, minPmaxT, ...)
})
################################################################################
# Coerce numeric classlabel to be integer, pass-on
#' @aliases mt,otuTable,numeric-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "numeric"), function(physeq, classlabel, minPmaxT="minP", ...){
	mt(physeq, as(classlabel, "integer"), minPmaxT="minP", ...)
})
################################################################################
# Coerce logical to integer, pass-on
#' @aliases mt,otuTable,logical-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "logical"), function(physeq, classlabel, minPmaxT="minP", ...){
	mt(physeq, as(classlabel, "integer"), minPmaxT="minP", ...)
})
################################################################################
# Test for length, then dispatch...
#' @aliases mt,otuTable,character-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "character"), function(physeq, classlabel, minPmaxT="minP", ...){
	if( length(classlabel) != nsamples(physeq) ){
		stop("classlabel not the same length as nsamples")
	} else {
		classlabel <- factor(classlabel)
	}
	# Use mt dispatch with classlabel now a suitable classlabel
	mt(physeq, classlabel, minPmaxT, ...)
})
################################################################################
#' @aliases mt,otuTable,factor-method
#' @rdname mt-methods
setMethod("mt", c("otuTable", "factor"), function(physeq, classlabel, minPmaxT="minP", ...){
	if( length(levels(classlabel)) != 2 ){
		stop("classlabel argument should be (or specify in sampleData) a two-level factor.\n",
		"In this case, the classlabel had ", length(levels(classlabel)), " levels (after coercion)."
		)
	}
	# integerize classlabel to logical suitable for mtmaxT
	classlabel <- (0:1)[classlabel]
	# Use mt dispatch with classlabel now a suitable classlabel
	mt(physeq, classlabel, minPmaxT, ...)
})
####################################################################################
# Internal function
# @aliases mt,matrix,integer-method
# not exported
# @keywords internal
mt.phyloseq.internal <- function(physeq, classlabel, minPmaxT="minP", ...){
	# require(multtest)
	if( minPmaxT == "minP" ){
		mt.minP(physeq, classlabel, ...)
	} else if( minPmaxT == "maxT" ){
		mt.maxT(physeq, classlabel, ...)	
	} else {
		print("Nothing calculated. minPmaxT argument must be either minP or maxT.")
	}
}
####################################################################################