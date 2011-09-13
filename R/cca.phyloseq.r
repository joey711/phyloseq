################################################################################
# vegan::cca "extension".
# formula is main input to this function. This complicates signature handling.
# A new method with a separate name is defined instead.
#
# Must transpose the phyloseq otuTable to fit the vegan::cca convention
# Whether-or-not to transpose needs to be a check, based on the 
#   "SpeciesAreRows" slot value
################################################################################
#' phyloseq package extension for \code{vegan::cca}.
#'
#' A formula is main input to \code{vegan::cca}. This complicates dispatch based
#' on object signature. A new method with a separate name is defined instead.
#' 
#' @param X formula-class object specifying the input. No need to splat higher-order
#'  objects. \code{cca.phyloseq} understands where to find the abundance table
#'  and sample characteristics.
#'
#' @param data \code{data.frame} containing information equivalent to
#'  a \code{sampleMap} object / component. Only necessary if complex object
#'  does not already contain \code{sampleMap} or you are keeping the data 
#'  separate for some reason.
#'
#' @return same output as \code{vegan::cca}.
#'
#' @seealso \code{\link{rda.phyloseq}}
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
####################################################################################
# Have not implemented cca.phyloseq for the (X, Y, Z, ...) variant of cca inputs. yet.
####################################################################################
# setMethod("cca.phyloseq", "otuTable", function(X, ...){
	# require(vegan)
	# if( speciesAreRows(X) ){ X <- t(X)@.Data }
	# vegan::cca(X=X, ...)
# })
# ####################################################################################
# setMethod("cca.phyloseq", "phyloseq", function(X, Y, ...){
	# # Option to specify which variables to contrain-by, rather than provide the whole matrix
	# if(class(Y) == "character"){ Y <- data.frame(sampleMap(X))[, Y, drop=FALSE] }
	# if(class(Z) == "character"){ Z <- data.frame(sampleMap(X))[, Z, drop=FALSE] }
	# cca.phyloseq(otuTable(X), Y, Z, ...)
# })
####################################################################################
# examples
####################################################################################
# cca.phyloseq(ex2 ~ Diet + Sex)
