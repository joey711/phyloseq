################################################################################
#' Build or access the taxonomyTable.
#'
#' This is the suggested method for both constructing and accessing a table of
#' taxonomic names, organized with ranks as columns (\code{\link{taxonomyTable-class}}). 
#' When the argument is a character matrix, tax_table() will create and return a 
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
#' @usage tax_table(object, errorIfNULL=TRUE)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain taxonomyTable.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return A \code{\link{taxonomyTable-class}} object.
#' It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is a 
#' character matrix representing the taxonomic classification of 
#' species in the experiment.
#'
#' @seealso \code{\link{phy_tree}}, \code{\link{sample_data}}, \code{\link{otu_table}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#'
#' @rdname tax_table-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # tax1 <- tax_table(matrix("abc", 30, 8))
#' # data(GlobalPatterns)
#' # tax_table(GlobalPatterns)
setGeneric("tax_table", function(object, errorIfNULL=TRUE) standardGeneric("tax_table"))
#' @rdname tax_table-methods
#' @aliases tax_table,ANY-method
setMethod("tax_table",  "ANY", function(object, errorIfNULL=TRUE){
	access(object, "tax_table", errorIfNULL)
})
# Constructor; for creating taxonomyTable from a matrix.
#' @rdname tax_table-methods
#' @aliases tax_table,matrix-method
setMethod("tax_table", "matrix", function(object){
  # Want dummy species/taxa index names if missing
  if(is.null(rownames(object))){
    rownames(object) <- paste("sp", 1:nrow(object), sep="")
  }
  if(is.null(colnames(object))){
    colnames(object) <- paste("ta", 1:ncol(object), sep="")
  }	
	# instantiate as taxonomyTable
	return(new("taxonomyTable", object))
})
# Constructor; coerce to matrix, then pass on for creating taxonomyTable.
#' @rdname tax_table-methods
#' @aliases tax_table,data.frame-method
setMethod("tax_table", "data.frame", function(object){
	# Warn first
  text = "Coercing from data.frame class to character matrix \n"
  text = paste0(text, "prior to building taxonomyTable. \n")
  text = paste0(text, "This could introduce artifacts. \n")
  text = paste0(text, "Check your taxonomyTable, or coerce to matrix manually.")
	warning(text)
	# Coerce everything to a matrix, then char-vector, then back to matrix.
	TT <- matrix(as(as(object, "matrix"), "character"),
               nrow=nrow(object),
               ncol=ncol(object)
        )
	# Pass on to matrix-method.
	tax_table(TT)
})
################################################################################
#' Subset species by taxonomic expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to speed subsetting complex experimental objects with one
#' function call. In the case of \code{subset_taxa}, the subsetting will be
#' based on an expression related to the columns and values within the 
#' \code{tax_table} (\code{taxonomyTable} component) slot of \code{physeq}.
#'
#' @usage subset_taxa(physeq, ...)
#'
#' @param physeq A \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}} that contains a
#'  taxonomyTable. If the \code{tax_table} slot is missing in \code{physeq}, then \code{physeq}
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
#' @rdname subset_taxa-methods
#' @docType methods
#' @export
#'
#' @examples
#' ## ex3 <- subset_taxa(GlobalPatterns, Phylum=="Bacteroidetes")
subset_taxa <- function(physeq, ...){
	if( is.null(tax_table(physeq)) ){ 
		cat("Nothing subset. No taxonomyTable in physeq.\n")
		return(physeq)
	} else {
		oldMA <- as(tax_table(physeq), "matrix")
		oldDF <- data.frame(oldMA)
		newDF <- subset(oldDF, ...)
		newMA <- as(newDF, "matrix")
		if( inherits(physeq, "taxonomyTable") ){
			return(tax_table(newMA))
		} else {
			tax_table(physeq) <- tax_table(newMA)
			return(physeq)
		}
	}
}
################################################################################
