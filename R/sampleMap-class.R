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