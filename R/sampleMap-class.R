################################################################################
#' Cleans absent levels in sampleMap/data.frame.
#'
#' This is used internally by the builder method, \code{\link{sampleMap}}, to
#' ensure that the factors describing categorical variables in a data.frame or
#' sampleMap object are free of extra levels that can plague downstream plots 
#' analysis.
#'
#' @usage reconcile_categories(DFSM)
#'
#' @param DFSM (Required). A \code{data.frame} or \code{sampleMap} object that needs to be cleaned. 
#'
#' @return A single \code{data.frame} object. Even if the input argument is a \code{sampleMap},
#'  the return is a \code{data.frame}. Because this is intended to be used internally by
#'  the builder method, it cannot also call the builder function to re-build
#'  the cleaned \code{sampleMap}.
#'
#' @keywords internal
#'
#' @examples
#' # # # data(ex1)
#' # # # SM <- sampleMap(ex1)
#' # # # DF <- data.frame(SM)
#' # # # DF <- data.frame(DF, col1=1:nrow(DF), col2=paste(1:nrow(DF), "t", sep=""))
#' # # # DF <- reconcile_variables(DF)
#' # # # SM <- sampleMap(reconcile_variables(SM))
#' # # # sapply(DF, class)
#' # # # sapply(SM, class)
reconcile_categories <- function(DFSM){
	DF <- as(DFSM, "data.frame")
	variable_classes <- sapply(DF, class)
	factor_cols <- names(variable_classes[variable_classes %in% c("factor", "character")])
	for( j in factor_cols){
		DF[, j] <- factor( as(DF[, j], "character") )
	}
	return( DF )
}
################################################################################
#' Subset samples by sampleMap expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to allow subsetting complex experimental objects with one
#' function call. In the case of \code{subset_samples}, the subsetting will be
#' based on an expression related to the columns and values within the 
#' sampleMap slot of \code{x}.
#'
#' @usage subset_samples(x, ...)
#'
#' @param x A \code{sampleMap}, or a \code{phyloseq-class} object with a
#'  \code{sampleMap}. If the \code{sampleMap} slot is missing in \code{x},
#'  then \code{x} will be returned as-is, and a warning will be printed to screen.
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
		cat("Nothing subset. No sampleMap in x.\n")
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