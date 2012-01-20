################################################################################
# Extend vegdist for phyloseq classes
# @importFrom vegan vegdist
################################################################################
#' \code{\link[vegan]{vegdist}} wrapper for phyloseq classes
#'
#' Trivially-extended S4 method from the \code{\link[vegan]{vegdist}} function,
#' such that S4 classes from the \code{\link{phyloseq-package}} are properly
#' handled / accessed. All parameters passed on to \code{\link[vegan]{vegdist}}
#' verbatim.
#'
#' @seealso \code{\link[vegan]{vegdist}} 
#' @import vegan
#' @export
#' @rdname vegdist-methods
#' @docType methods
#' @aliases vegdist
#'
#' @examples
#' # data(esophagus)
#' # vegdist(esophagus, "jaccard")
setGeneric("vegdist")
################################################################################
#' @aliases vegdist,otuTable-method
#' @rdname vegdist-methods
setMethod("vegdist", "otuTable", function(x, method = "bray", binary = FALSE,
	diag = FALSE, upper = FALSE, na.rm = FALSE, ...){
	# Make sure in sample-by-species orientation
	if( speciesAreRows(x) ){x <- t(x)}
	# Convert to simple matrix
	x <- as(x, "matrix")
	# pass to standard method (compiled C)
	vegdist(x, method, binary, diag, upper, na.rm, ...)	
})
################################################################################
#' @aliases vegdist,phyloseq-method
#' @rdname vegdist-methods
setMethod("vegdist", "phyloseq", function(x, method = "bray", binary = FALSE,
	diag = FALSE, upper = FALSE, na.rm = FALSE, ...){
	# Simply access the otuTable
	x <- otuTable(x)
	vegdist(x, method, binary, diag, upper, na.rm, ...)	
})
################################################################################