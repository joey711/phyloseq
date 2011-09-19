################################################################################
#' Subset samples by sampleMap expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to allow subsetting complex experimental objects with one
#' function call. In the case of \code{subset_samples}, the subsetting will be
#' based on an expression related to the columns and values within the 
#' sampleMap slot of \code{x}.
#'
#' @param x A \code{sampleMap} class, or a more complex class that contains a
#'  \code{sampleMap}. If the \code{sampleMap} slot is missing in \code{x}, then \code{x}
#'  will be returned as-is and a warning will be printed to screen.
#'
#' @param ... The subsetting expression that should be applied to the 
#'  \code{sampleMap}. This is passed on to \code{\link{subset}}, see its
#'  documentation for more details.
#'
#' @return A subsetted object with the same class as \code{x}.
#' 
#' @seealso \code{\link{subset_species}}
#'
#' @export
#'
#' @examples
#' ## data(ex1)
#' ## ex2 <- subset_samples(ex1, Gender=="A")
#' ## ex2 <- subset_samples(ex1, Diet==1)
subset_samples <- function(x, ...){
	if( is.null(sampleMap(x)) ){ 
		cat("Nothing subset. No TaxonomyTable in x.\n")
		return(x)
	} else {
		oldDF <- as(sampleMap(x), "data.frame")
		newDF <- subset(oldDF, ...)
		if( class(x) == "sampleMap" ){
			return(sampleMap(newDF))
		} else {
			sampleMap(x) <- sampleMap(newDF)
			return(x)
		}
	}
}
################################################################################
#' Subset species by taxonomic expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to speed subsetting complex experimental objects with one
#' function call. In the case of \code{subset_species}, the subsetting will be
#' based on an expression related to the columns and values within the 
#' taxTab slot of \code{x}.
#'
#' @param x A taxonomyTable class, or a more complex class that contains a
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
		cat("Nothing subset. No TaxonomyTable in x.\n")
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