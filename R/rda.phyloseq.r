####################################################################################
#' phyloseq package extension for \code{vegan::rda}.
#'
#' A formula is main input to \code{vegan::cca}. This complicates dispatch based
#' on object signature. A new method with a separate name is defined instead.
#' 
#' @param X formula-class object specifying the input. No need to splat higher-order
#'  objects. \code{rda.phyloseq} understands where to find the abundance table
#'  and sample characteristics.
#'
#' @param data \code{data.frame} containing information equivalent to
#'  a \code{sampleMap} object / component. Only necessary if complex object
#'  does not already contain \code{sampleMap} or you are keeping the data 
#'  separate for some reason.
#'
#' @return same output as \code{\link{rda}}
#' @keywords redundancy analysis RDA rda
#' @seealso cca.phyloseq
#' @import vegan
#' @export
#' @examples #
#' ## data(ex1)
#' ## rda.phyloseq(ex1 ~ Diet + Gender)
setGeneric("rda.phyloseq", function(X, ...) standardGeneric("rda.phyloseq"))
setMethod("rda.phyloseq", "formula", function(X, data=NULL){
	require(vegan)
	phylobject <- get( as.character(X)[2] )
	OTU        <- otuTable( phylobject )
	if( speciesAreRows(OTU) ){
		OTU <- t(OTU)@.Data
	} else {
		OTU <- OTU@.Data
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sampleMap slot
	if( is.null(data) ){
		data <- data.frame(sampleMap(phylobject))
	}
	vegan::rda(newFormula, data=data)	
})
####################################################################################