####################################################################################
# # # # Avoiding full import of multtest to mitigate potential conflicts
####################################################################################
#' Multiple testing of taxa abundance acccording to sample categories/classes
#'
#' @usage mt(physeq, classlabel, minPmaxT="minP", ...)
#'
#' @param physeq (Required). \code{\link{otu_table-class}} or \code{\link{phyloseq-class}}.
#'  In this multiple testing framework, different taxa correspond to variables
#'  (hypotheses), and samples to observations.
#'
#' @param classlabel (Required). A single character index of the sample-variable
#'  in the \code{\link{sample_data}} of \code{physeq} that will be used for multiple testing. 
#'  Alternatively, \code{classlabel} can be a custom integer (or numeric coercable
#'  to an integer), character, or factor with
#'  length equal to \code{nsamples(physeq)}. 
#'  
#'  NOTE: the default test applied to each taxa is a two-sample two-sided 
#'  \code{\link{t.test}}, WHICH WILL FAIL with an error if you provide a data variable
#'  (or custom vector) that contains MORE THAN TWO classes. One alternative to consider
#'  is an F-test, by specifying \code{test="f"} as an additional argument. See
#'  the first example below, and/or further documentation of
#'  \code{\link[multtest]{mt.maxT}} or \code{\link[multtest]{mt.minP}}
#'  for other options and formal details.
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
#' ## # Simple example, testing genera that sig correlate with Enterotypes
#' ## data(enterotype)
#' ## # Filter samples that don't have Enterotype
#' ## x <- subset_samples(enterotype, !is.na(Enterotype))
#' ## # (the taxa are at the genera level in this dataset)
#' ## mt(x, "Enterotype", test="f")
#' ## # Not surprisingly, Prevotella and Bacteroides top the list.
#' ## # Different test, multiple-adjusted t-test, whether samples are ent-2 or not.
#' ## mt(x, get_variable(x, "Enterotype")==2)
setGeneric("mt", function(physeq, classlabel, minPmaxT="minP", ...) standardGeneric("mt") )
################################################################################
# First, access the otu_table, and if appropriate, define classlabel from 
# the sample_data.
#' @aliases mt,phyloseq,ANY-method
#' @rdname mt-methods
setMethod("mt", c("phyloseq", "ANY"), function(physeq, classlabel, minPmaxT="minP", ...){
	# Extract the class information from the sample_data
	# if sample_data slot is non-empty,
	# and the classlabel is a character-class
	# and its length is 1.
	if( !is.null(sample_data(physeq, FALSE)) & 
			inherits(classlabel, "character") &
			identical(length(classlabel), 1L) ){
		# Define a raw factor based on the data available in a sample variable
		rawFactor = get_variable(physeq, classlabel[1])
		if( !inherits(rawFactor, "factor") ){
			# coerce to a factor if it is not already one.
			rawFactor = factor(rawFactor)
		}
		# Either way, replace `classlabel` with `rawFactor`
		classlabel = rawFactor
	}
	# Either way, dispatch `mt` on otu_table(physeq)
	MT = mt(otu_table(physeq), classlabel, minPmaxT, ...)
	if( !is.null(tax_table(physeq, FALSE)) ){
		# If there is tax_table data present,
		# add/cbind it to the results.		
		MT = cbind(MT, as(tax_table(physeq), "matrix")[rownames(MT), ])
	}
	return(MT)
})
################################################################################
# All valid mt() calls eventually funnel dispatch to this method.
# The otu_table orientation is checked/handled here (and only here).
#' @aliases mt,otu_table,integer-method
#' @rdname mt-methods
setMethod("mt", c("otu_table", "integer"), function(physeq, classlabel, minPmaxT="minP", ...){
	# Guarantee proper orientation of abundance table, and coerce to matrix.
	if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)	}
	mt.phyloseq.internal(as(physeq, "matrix"), classlabel, minPmaxT, ...)
})
################################################################################
# Coerce numeric classlabel to be integer, pass-on
#' @aliases mt,otu_table,numeric-method
#' @rdname mt-methods
setMethod("mt", c("otu_table", "numeric"), function(physeq, classlabel, minPmaxT="minP", ...){
	mt(physeq, as(classlabel, "integer"), minPmaxT="minP", ...)
})
################################################################################
# Coerce logical to integer, pass-on
#' @aliases mt,otu_table,logical-method
#' @rdname mt-methods
setMethod("mt", c("otu_table", "logical"), function(physeq, classlabel, minPmaxT="minP", ...){
	mt(physeq, as(classlabel, "integer"), minPmaxT="minP", ...)
})
################################################################################
# Test for length, then dispatch...
#' @aliases mt,otu_table,character-method
#' @rdname mt-methods
setMethod("mt", c("otu_table", "character"), function(physeq, classlabel, minPmaxT="minP", ...){
	if( length(classlabel) != nsamples(physeq) ){
		stop("classlabel not the same length as nsamples(physeq)")
	} else {
		classlabel <- factor(classlabel)
	}
	# Use mt dispatch with classlabel now a suitable classlabel
	mt(physeq, classlabel, minPmaxT, ...)
})
################################################################################
# Coerce factor to an integer vector of group labels,
# starting at 0 for the first group
#' @aliases mt,otu_table,factor-method
#' @rdname mt-methods
setMethod("mt", c("otu_table", "factor"), function(physeq, classlabel, minPmaxT="minP", ...){
	# integerize classlabel, starting at 0
	classlabel <- (0:(length(classlabel)-1))[classlabel]
	# Use mt dispatch with classlabel now a suitable classlabel
	mt(physeq, classlabel, minPmaxT, ...)
})
####################################################################################
# Internal function
# @aliases mt,matrix,integer-method
# not exported
#' @keywords internal
mt.phyloseq.internal <- function(physeq, classlabel, minPmaxT="minP", ...){
	# require(multtest)
	if( minPmaxT == "minP" ){
		return( mt.minP(physeq, classlabel, ...) )
	} else if( minPmaxT == "maxT" ){
		return( mt.maxT(physeq, classlabel, ...) )
	} else {
		print("Nothing calculated. minPmaxT argument must be either minP or maxT.")
	}
}
####################################################################################
