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