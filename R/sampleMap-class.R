################################################################################
#' Build or access the sampleMap.
#'
#' This is the suggested method for both constructing and accessing a table
#' of sample-level variables (\code{\link{sampleMap-class}}), 
#' which in the \code{\link{phyloseq-package}} is represented as a special
#' extension of the \code{\link{data.frame-class}}.
#' When the
#' argument is a data.frame, sampleMap() will create a sampleMap-class object.
#' In this case, the rows should be named to match the
#' \code{sample.names} of the other objects to which it will ultimately be paired.
#' Alternatively, if the first argument is an experiment-level (\code{\link{phyloseq-class}})
#' object, then the corresponding \code{sampleMap} is returned.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{phyloseq-class} but does not 
#' contain a \code{sampleMap}.
#'
#' @usage sampleMap(object, errorIfNULL=TRUE)
#'
#' @param object (Required). A \code{\link{data.frame-class}}, 
#'  or a \code{\link{phyloseq-class}} object.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. 
#'
#' @return A \code{\link{sampleMap-class}} object
#' representing the sample variates of an experiment.
#'
#' @seealso \code{\link{tre}}, \code{\link{taxTab}}, \code{\link{otuTable}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#'
#' @aliases sampleMap samplemap
#'
#' @rdname sampleMap-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # data(ex1)
#' # sampleMap(ex1)
setGeneric("sampleMap", function(object, errorIfNULL=TRUE) standardGeneric("sampleMap"))
#' @rdname sampleMap-methods
#' @aliases sampleMap,ANY-method
setMethod("sampleMap", "ANY", function(object, errorIfNULL=TRUE){
	access(object, "sampleMap", errorIfNULL)
})
# constructor; for creating sampleMap from a data.frame
#' @rdname sampleMap-methods
#' @aliases sampleMap,data.frame-method
setMethod("sampleMap", "data.frame", function(object){
	# Make sure there are no phantom levels in categorical variables
	object <- reconcile_categories(object)

	# instantiate first to check validity
	SM <- new("sampleMap", object)
		
	# Want dummy samples index names if missing
	if( all(rownames(SM) == as.character(1:nrow(SM))) ){
		rownames(SM) <- paste("sa", 1:nrow(SM), sep="")
	}	
	return(SM)
})
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
#' function call. The subsetting will be
#' based on an expression related to the columns and values within the 
#' sampleMap.
#'
#' @usage subset_samples(physeq, ...)
#'
#' @param physeq A \code{\link{sampleMap-class}}, or a \code{\link{phyloseq-class}}
#'  object with a 
#'  \code{sampleMap}. If the \code{sampleMap} slot is missing in \code{physeq},
#'  then \code{physeq} will be returned as-is, and a warning will be printed to screen.
#'
#' @param ... The subsetting expression that should be applied to the 
#'  \code{sampleMap}. This is passed on to \code{\link{subset}}, see its
#'  documentation for more details.
#'
#' @return A subsetted object with the same class as \code{physeq}.
#' 
#' @seealso \code{\link{subset_species}}
#'
#' @export
#' @rdname subset_samples-methods
#' @docType methods
#'
#' @examples
#'  # data(ex1)
#'  # ex2 <- subset_samples(ex1, Gender=="A")
#'  # ex2 <- subset_samples(ex1, Diet==1)
#'  ### Here is an example comparing subset_samples with prune_samples...
#'  # B_only_sample_names <- sample.names(sampleMap(ex1)[(sampleMap(ex1)[, "Gender"]=="B"),])
#'  # ex2 <- prune_samples(B_only_sample_names, ex1)
#'  # ex3 <- subset_samples(ex1, Gender=="B")
#'  # ## This should be TRUE.
#'  # identical(ex2, ex3)
subset_samples <- function(physeq, ...){
	if( is.null(sampleMap(physeq)) ){ 
		cat("Nothing subset. No sampleMap in physeq.\n")
		return(physeq)
	} else {
		oldDF <- as(sampleMap(physeq), "data.frame")
		newDF <- subset(oldDF, ...)
		if( class(physeq) == "sampleMap" ){
			return(sampleMap(newDF))
		} else {
			sampleMap(physeq) <- sampleMap(newDF)
			return(physeq)
		}
	}
}
################################################################################