################################################################################
# vegan::cca "extension".
# formula is main input to this function. This complicates signature handling.
# A new method with a separate name is defined instead.
#
# Must transpose the phyloseq otuTable to fit the vegan::cca convention
# Whether-or-not to transpose needs to be a check, based on the 
#   "SpeciesAreRows" slot value
################################################################################
#' phyloseq package extension for \code{\link[vegan]{cca}}.
#'
#' A formula is main input to \code{\link[vegan]{cca}}. This complicates dispatch based
#' on object signature. A new method with a separate name is defined instead.
#'
#' @usage cca.phyloseq(X, ...)
#' 
#' @param X (Required). A formula-class object, specifying the input.
#'  No need to splat higher-order
#'  objects. \code{cca.phyloseq} understands where to find the abundance table
#'  and sample characteristics. Alternatively, \code{X} can be an otuTable, or any
#'  more complex phyloseq-package class that contains an otuTable. 
#'
#' @param ... (Optional). E.g. \code{data=DF}, where \code{DF} is a \code{data.frame}
#'  containing information equivalent to
#'  a \code{sampleMap} object / component. Only necessary if complex object
#'  does not already contain \code{sampleMap} or you are keeping the data 
#'  separate for some reason.
#'
#' @return same output as \code{\link[vegan]{cca}}.
#'
#' @seealso \code{\link{rda.phyloseq}}, \code{\link[vegan]{cca}}
#'
#' @rdname cca.phyloseq-methods
#' @docType methods
#'
#' @export
#' @import vegan
#' @examples #
#' ## data(ex1)
#' ## cca.phyloseq(ex1 ~ Diet + Gender)
setGeneric("cca.phyloseq", function(X, ...) standardGeneric("cca.phyloseq"))
################################################################################
#' @aliases cca.phyloseq,formula-method
#' @rdname cca.phyloseq-methods
setMethod("cca.phyloseq", "formula", function(X, data=NULL){
	phylobject <- get( as.character(X)[2] )
	OTU        <- otuTable( phylobject )
	if( speciesAreRows(OTU) ){
		OTU <- t(as(OTU, "matrix"))
	} else {
		OTU <- as(OTU, "matrix")
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sampleMap slot
	if( is.null(data) ){
		data <- data.frame(sampleMap(phylobject))
	}
	vegan::cca(newFormula, data=data)	
})
################################################################################
#' @aliases cca.phyloseq,otuTable-method
#' @rdname cca.phyloseq-methods
setMethod("cca.phyloseq", "otuTable", function(X){
	if( speciesAreRows(X) ){
		X <- t(as(X, "matrix"))
	} else {
		X <- as(X, "matrix")
	}
	vegan::cca(X)	
})
################################################################################
#' @aliases cca.phyloseq,phyloseq-method
#' @rdname cca.phyloseq-methods
setMethod("cca.phyloseq", "phyloseq", function(X){
	cca.phyloseq(otuTable(X))
})
################################################################################
