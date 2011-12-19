################################################################################
#' Build or access the taxonomyTable.
#'
#' This is the suggested method for both constructing and accessing a table of
#' taxonomic names, organized with ranks as columns (\code{\link{taxonomyTable-class}}). 
#' When the argument is a character matrix, taxTab() will create and return a 
#' \code{\link{taxonomyTable-class}} object.
#' In this case, the rows should be named to match the
#' \code{species.names} of the other objects to which it will ultimately be paired.
#' Alternatively, if the first argument is an experiment-level (\code{\link{phyloseq-class}})
#' object, then the corresponding \code{taxonomyTable} is returned.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{phyloseq-class} but does not 
#' contain a \code{taxonomyTable}. 
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
#' @seealso \code{\link{tre}}, \code{\link{sampleData}}, \code{\link{otuTable}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#'
#' @rdname taxTab-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # data(ex1)
#' # taxTab(ex1)
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
	# instantiate first to check validity
	TT <- new("taxonomyTable", object)
		
	# Want dummy species/taxa index names if missing
	if(is.null(rownames(TT))){
		rownames(TT) <- paste("sp", 1:nrow(TT), sep="")
	}
	if(is.null(colnames(TT))){
		colnames(TT) <- paste("ta", 1:ncol(TT), sep="")
	}	
	return(TT)
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
#' \code{taxTab} (\code{taxonomyTable} component) slot of \code{physeq}.
#'
#' @usage subset_species(physeq, ...)
#'
#' @param physeq A \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}} that contains a
#'  taxonomyTable. If the \code{taxTab} slot is missing in \code{physeq}, then \code{physeq}
#'  will be returned as-is and a warning will be printed to screen.
#'
#' @param ... The subsetting expression that should be applied to the 
#'  \code{taxonomyTable}. This is passed on to \code{\link{subset}}, and more
#'  details and examples about how it functions can be found in its documentation.
#'
#' @return A subsetted object with the same class as \code{physeq}.
#' 
#' @seealso \code{\link{subset_samples}}
#'
#' @rdname subset_species-methods
#' @docType methods
#' @export
#'
#' @examples
#' ## ex3 <- subset_species(ex1, Phylum=="Bacteroidetes")
subset_species <- function(physeq, ...){
	if( is.null(taxTab(physeq)) ){ 
		cat("Nothing subset. No taxonomyTable in physeq.\n")
		return(physeq)
	} else {
		oldMA <- as(taxTab(physeq), "matrix")
		oldDF <- data.frame(oldMA)
		newDF <- subset(oldDF, ...)
		newMA <- as(newDF, "matrix")
		if( class(physeq) == "taxonomyTable" ){
			return(taxTab(newMA))
		} else {
			taxTab(physeq) <- taxTab(newMA)
			return(physeq)
		}
	}
}
################################################################################