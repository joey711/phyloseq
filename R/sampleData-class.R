################################################################################
#' Build or access sample_data.
#'
#' This is the suggested method for both constructing and accessing a table
#' of sample-level variables (\code{\link{sample_data-class}}), 
#' which in the \code{\link{phyloseq-package}} is represented as a special
#' extension of the \code{\link{data.frame-class}}.
#' When the
#' argument is a \code{\link{data.frame}}, \code{sample_data} will create 
#' a sample_data-class object.
#' In this case, the rows should be named to match the
#' \code{\link{sample_names}} of the other objects to which it will ultimately be paired.
#' Alternatively, if the first argument is an experiment-level (\code{\link{phyloseq-class}})
#' object, then the corresponding \code{sample_data} is returned.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{phyloseq-class} but does not 
#' contain a \code{sample_data}.
#'
#' @usage sample_data(object, errorIfNULL=TRUE)
#'
#' @param object (Required). A \code{\link{data.frame-class}}, 
#'  or a \code{\link{phyloseq-class}} object.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. 
#'
#' @return A \code{\link{sample_data-class}} object
#' representing the sample variates of an experiment.
#'
#' @seealso \code{\link{phy_tree}}, \code{\link{tax_table}}, \code{\link{otu_table}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#'
#' @aliases sample_data
#'
#' @rdname sample_data-methods
#' @docType methods
#' @export
#'
#' @examples #
#' data(soilrep)
#' head(sample_data(soilrep))
setGeneric("sample_data", function(object, errorIfNULL=TRUE) standardGeneric("sample_data"))
#' @rdname sample_data-methods
#' @aliases sample_data,ANY-method
setMethod("sample_data", "ANY", function(object, errorIfNULL=TRUE){
	access(object, "sam_data", errorIfNULL)
})
# constructor; for creating sample_data from a data.frame
#' @rdname sample_data-methods
#' @aliases sample_data,data.frame-method
setMethod("sample_data", "data.frame", function(object){
	# Make sure there are no phantom levels in categorical variables
	object <- reconcile_categories(object)

	# instantiate first to check validity
	SM <- new("sample_data", object)
		
	# Want dummy samples index names if missing
	if( all(rownames(SM) == as.character(1:nrow(SM))) ){
		rownames(SM) <- paste("sa", 1:nrow(SM), sep="")
	}	
	return(SM)
})
################################################################################
#' Cleans absent levels in sample_data/data.frame.
#'
#' This is used internally by the builder method, \code{\link{sample_data}}, to
#' ensure that the factors describing categorical variables in a data.frame or
#' sample_data object are free of extra levels that can plague downstream plots 
#' analysis.
#'
#' @usage reconcile_categories(DFSM)
#'
#' @param DFSM (Required). A \code{data.frame} or \code{sample_data} object that needs to be cleaned. 
#'
#' @return A single \code{data.frame} object. Even if the input argument is a \code{sample_data},
#'  the return is a \code{data.frame}. Because this is intended to be used internally by
#'  the builder method, it cannot also call the builder function to re-build
#'  the cleaned \code{sample_data}.
#'
#' @keywords internal
#'
#' @examples
#' # # # data(GlobalPatterns)
#' # # # SM <- sample_data(GlobalPatterns)
#' # # # DF <- data.frame(SM)
#' # # # DF <- data.frame(DF, col1=1:nrow(DF), col2=paste(1:nrow(DF), "t", sep=""))
#' # # # DF <- reconcile_categories(DF)
#' # # # SM <- sample_data(reconcile_categories(SM))
#' # # # sapply(DF, class)
#' # # # sapply(SM, class)
reconcile_categories <- function(DFSM){
	DF = as(DFSM, "data.frame")
	#variable_classes <- sapply(DF, class)
	#factor_cols <- names(variable_classes[variable_classes %in% c("factor", "character")])
	factor_cols = which(sapply(DF, inherits, what="factor"))
	for( j in factor_cols){
		DF[, j] <- factor( as(DF[, j], "character") )
	}
	return(DF)
}
################################################################################
#' Subset samples by sample_data expression
#'
#' This is a convenience wrapper around the \code{\link{subset}} function.
#' It is intended to allow subsetting complex experimental objects with one
#' function call.
#' Subsetting is based on an expression for which the context first includes
#' the variables contained in \code{\link{sample_data}}.
#' The \code{samples} retained in the dataset is equivalent to
#' \code{x[subset & !is.na(subset)]}, where \code{x} is the vector of sample IDs
#' and \code{subset} is the logical that results from your subsetting expression.
#' This is important to keep in mind, as users are often unaware that this
#' subsetting step also removes/omits samples that have a missing value, \code{NA},
#' somewhere in the expression.
#'
#' @usage subset_samples(physeq, ...)
#'
#' @param physeq A \code{\link{sample_data-class}}, or a \code{\link{phyloseq-class}}
#'  object with a 
#'  \code{sample_data}. If the \code{sample_data} slot is missing in \code{physeq},
#'  then \code{physeq} will be returned as-is, and a warning will be printed to screen.
#'
#' @param ... The subsetting expression that should be applied to the 
#'  \code{sample_data}. This is passed on to \code{\link{subset}}, see its
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
#'  # data(GlobalPatterns)
#'  # subset_samples(GlobalPatterns, SampleType=="Ocean")
subset_samples <- function(physeq, ...){
	if( is.null(sample_data(physeq)) ){ 
		cat("Nothing subset. No sample_data in physeq.\n")
		return(physeq)
	} else {
		oldDF <- as(sample_data(physeq), "data.frame")
		newDF <- subset(oldDF, ...)
		if( class(physeq) == "sample_data" ){
			return(sample_data(newDF))
		} else {
			sample_data(physeq) <- sample_data(newDF)
			return(physeq)
		}
	}
}
################################################################################
