################################################################################
#' Access taxTab slot, or instantiate taxonomyTable-class.
#'
#' \code{taxTab()} is both a constructor and accessor method. When the
#' argument is a character matrix, taxTab() will attempt to create and return a 
#' taxonomyTable-class object. In this case, the rows should be named to match the
#' \code{species.names} of the other objects to which it will ultimately be paired.
#' Alternatively, if the argument is an object that 
#' contains a taxonomyTable, including a taxonomyTable-class object, then the 
#' corresponding taxonomyTable is returned, as-is.
#'
#' This is the main method suggested for constructing taxonomyTable objects from 
#' a character matrix of taxonomy classifications. Each row represents a 
#' different species.
#' It is also the suggested method for accessing/subsetting
#' the taxonomyTable from a more complex object. \code{taxTab} is the slot name
#' that holds the taxonomyTable-class object in a multi-component phyloseq
#' object.
#'
#' @usage taxTab(object, errorIfNULL=TRUE)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain taxonomyTable.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return A taxonomyTable object. It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is a 
#' character matrix representing the taxonomic classification of species in the
#' experiment.
#'
#' @seealso otuTable sampleMap tre phyloseq
#' @aliases taxTab taxtab
#'
#' @rdname taxTab-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # taxTab(ex1)
#' # tax1
setGeneric("taxTab", function(object, errorIfNULL=TRUE) standardGeneric("taxTab"))
#' @rdname taxTab-methods
#' @aliases taxTab,ANY-method
setMethod("taxTab",  "ANY", function(object, errorIfNULL=TRUE){
	access(object, "taxTab", errorIfNULL)
})
# Constructor; for creating taxonomyTable from a matrix.
#' @rdname taxTab-methods
#' @aliases taxTab,matrix-method
setMethod("taxTab", "matrix", function(object){
	# Want dummy species/taxa index names if missing
	if(is.null(rownames(object))){
		rownames(object) <- paste("sp", 1:nrow(object), sep="")
	}
	if(is.null(colnames(object))){
		colnames(object) <- paste("ta", 1:ncol(object), sep="")
	}	
	new("taxonomyTable", object)
})
#' @rdname taxTab-methods
#' @aliases taxTab taxtab
#' @export
taxtab <- taxTab
################################################################################
#' Subset species by taxonomic expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to speed subsetting complex experimental objects with one
#' function call. In the case of \code{subset_species}, the subsetting will be
#' based on an expression related to the columns and values within the 
#' \code{taxTab} (\code{taxonomyTable} component) slot of \code{x}.
#'
#' @usage subset_species(x, ...)
#'
#' @param x A taxonomyTable class, or \code{phyloseq-class} that contains a
#'  taxonomyTable. If the \code{taxTab} slot is missing in \code{x}, then \code{x}
#'  will be returned as-is and a warning will be printed to screen.
#'
#' @param ... The subsetting expression that should be applied to the 
#'  \code{taxonomyTable}. This is passed on to \code{\link{subset}}, and more
#'  details and examples about how it functions can be found in its documentation.
#'
#' @return A subsetted object with the same class as \code{x}.
#' 
#' @seealso \code{\link{subset_samples}}
#'
#' @export
#'
#' @examples
#' ## ex3 <- subset_species(ex1, Phylum=="Bacteroidetes")
subset_species <- function(x, ...){
	if( is.null(taxTab(x)) ){ 
		cat("Nothing subset. No taxonomyTable in x.\n")
		return(x)
	} else {
		oldMA <- as(taxTab(x), "matrix")
		oldDF <- data.frame(oldMA)
		newDF <- subset(oldDF, ...)
		newMA <- as(newDF, "matrix")
		if( class(x) == "taxonomyTable" ){
			return(taxTab(newMA))
		} else {
			taxTab(x) <- taxTab(newMA)
			return(x)
		}
	}
}
################################################################################