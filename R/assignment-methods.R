################################################################################
#' Assign a new OTU Table to \code{x}
#'
#' @usage otuTable(x) <- value
#'
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{otuTable-class}} or \code{\link{phyloseq-class}}.
#'
#' @export
#' @docType methods
#' @rdname assign-otuTable
#' @aliases assign-otuTable otuTable<-
#'
#' @examples
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # OTU <- otuTable(GlobalPatterns)[1:100, ]
#' # otuTable(ex2b) <- OTU
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Relace otuTable by implying the component in context.
#' # ex2c <- GlobalPatterns
#' # otuTable(ex2c) <- ex2b
#' # identical(ex2a, ex2c)
setGeneric("otuTable<-", function(x, value) standardGeneric("otuTable<-"))
#' @rdname assign-otuTable
#' @aliases otuTable<-,phyloseq,otuTable-method
setMethod("otuTable<-", c("phyloseq", "otuTable"), function(x, value){
	phyloseq(otuTable=value, samData=x@samData, taxTab=x@taxTab, tre=x@tre)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuTable,otuTable-method
setMethod("otuTable<-", c("otuTable", "otuTable"), function(x, value){ value })
#' @rdname assign-otuTable
#' @aliases otuTable<-,phyloseq,phyloseq-method
setMethod("otuTable<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(otuTable=otuTable(value), samData=x@samData, taxTab=x@taxTab, tre=x@tre)
})
################################################################################
#' Manually change speciesAreRows through assignment.
#'
#' The speciesAreRows slot is a logical indicating the orientation of the
#' abundance table contained in object \code{x}.
#'
#' @usage speciesarerows(x) <- value
#'
#' @param x \code{\link{otuTable-class}} or \code{\link{phyloseq-class}}
#'
#' @param value A logical of length equal to 1. If \code{length(value) > 1}, 
#'  the additional elements will be ignored. Only the first element is assigned
#'  to the speciesAreRows slot.
#'
#' @export
#' @docType methods
#' @rdname assign-speciesarerows
#' @aliases assign-speciesarerows speciesarerows<-
#'
#' @examples #
#' # data(GlobalPatterns)
#' # speciesarerows(GlobalPatterns)
#' # speciesarerows(otuTable(GlobalPatterns))
setGeneric("speciesarerows<-", function(x, value){
	standardGeneric("speciesarerows<-")
})
#' @rdname assign-speciesarerows
#' @aliases speciesarerows<-,otuTable,logical-method
setMethod("speciesarerows<-", c("otuTable", "logical"), function(x, value){
	x@speciesAreRows <- value[1]
	return(x)
})
#' @rdname assign-speciesarerows
#' @aliases speciesarerows<-,phyloseq,logical-method
setMethod("speciesarerows<-", c("phyloseq", "logical"), function(x, value){
	speciesarerows(otuTable(x)) <- value
	return(x)
})
################################################################################
#' Assign (new) sampleData to \code{x}
#'
#' This replaces the current \code{sampleData} component of \code{x} with 
#' \code{value}, if \code{value} is a \code{\link{sampleData-class}}. However,
#' if \code{value} is a \code{data.frame}, then \code{value} is first coerced to
#' a \code{\link{sampleData-class}}, and then assigned. Alternatively, if 
#' \code{value} is \code{\link{phyloseq-class}}, then the 
#' \code{\link{sampleData}} component will first be accessed from \code{value}
#'  and then assigned. This makes possible some concise assignment/replacement
#'  statements when adjusting, modifying, or building subsets of 
#'  experiment-level data. See some examples below.
#'  
#' Internally, this re-builds the \code{\link{phyloseq-class}} object using
#' the standard \code{\link{phyloseq}} constructor. Thus, index mismatches
#' between sample-describing components will not be allowed, and subsetting
#' will occurr automatically such that only the intersection of sample IDs
#' are included in any components. This has the added benefit of re-checking
#' (internally) for any other issues.
#'
#' @usage sampleData(x) <- value
#'
#' @param x (Required). \code{\link{phyloseq-class}}. The object to modify.
#' @param value (Required). Either a \code{\link{sampleData-class}}, 
#'  a \code{data.frame} that can be coerced into \code{\link{sampleData-class}}, 
#'  or a \code{\link{phyloseq-class}} that contains a 
#'  suitable \code{sampleData} component to assign to \code{x}. If unsure,
#'  try \code{\link{sampleData}}\code{(value)}, which should return a 
#'  \code{\link{sampleData-class}} object without error.
#'
#' @return No return. This is an assignment statement.
#'
#' @export
#' @rdname assign-sampleData
#' @aliases assign-sampleData sampleData<- samData<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 10 samples in GlobalPatterns
#' # ex2a <- prune_samples(sample.names(GlobalPatterns)[1:10], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # SD <- sampleData(GlobalPatterns)[1:10, ]
#' # sampleData(ex2b) <- SD
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Example restoring the original sampleData component. ex2c lacks sampleData
#' # ex2c <- phyloseq(otuTable(GlobalPatterns), taxTab(GlobalPatterns), tre(GlobalPatterns))
#' # sampleData(ex2c) <- GlobalPatterns
#' # identical(ex2c, GlobalPatterns)
#' # # Can try on ex2b, but other components have only 10 samples. No change.
#' # sampleData(ex2b) <- GlobalPatterns
#' # identical(ex2a, ex2b) # still true.
"sampleData<-" <- function(x, value){
	if(class(value) != "sampleData"){value <- sampleData(value)}
	phyloseq(otuTable=x@otuTable, samData=value, taxTab=x@taxTab, tre=x@tre)
}
#' @export
#' @rdname assign-sampleData
#' @aliases assign-sampleData sampleData<- samData<-
#' @usage samData(x) <- value
"samData<-" <- function(x, value){
	if(class(value) != "sampleData"){value <- sampleData(value)}
	phyloseq(otuTable=x@otuTable, samData=value, taxTab=x@taxTab, tre=x@tre)
}
################################################################################
#' Assign a (new) Taxonomy Table to \code{x}
#'
#' @usage taxTab(x) <- value
#' 
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{taxonomyTable-class}}.
#'  Alternatively, \code{value} can be a \code{\link{phyloseq-class}} that has
#'  a \code{\link{taxTab}} component, or a \code{\link{matrix-class}}
#'  that can be coerced to a \code{\link{taxonomyTable-class}} with row indices
#'  that match at least some of the \code{\link{species.names}} of \code{x}.
#'
#' @export
#' @rdname assign-taxTab
#' @aliases assign-taxTab taxTab<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # TT <- taxTab(GlobalPatterns)[1:100, ]
#' # taxTab(ex2b) <- TT
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # 2 examples adding a taxTab component from phyloseq or matrix classes
#' # ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), tre(ex2b))
#' # taxTab(ex2c) <- ex2b
#' # identical(ex2a, ex2c)
#' # ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), tre(ex2b))
#' # taxTab(ex2c) <- as(taxTab(ex2b), "matrix")
#' # identical(ex2a, ex2c)
"taxTab<-" <- function(x, value){
	if(class(value) != "taxonomyTable"){value <- taxTab(value)}
	phyloseq(otuTable=x@otuTable, samData=x@samData, taxTab=value, tre=x@tre)
}
################################################################################
#' Assign a (new) phylogenetic tree to \code{x}
#'
#' @usage tre(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{phylo-class}}, or \code{\link{phyloseq-class}}
#'
#' @export
#' @docType methods
#' @rdname assign-tre
#' @aliases assign-tre tre<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # tree <- prune_species(species.names(GlobalPatterns)[1:100], tre(GlobalPatterns))
#' # tre(ex2b) <- tree
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Example adding a phylo tree from phyloseq class
#' # ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), taxTab(ex2b))
#' # tre(ex2c) <- ex2b
#' # identical(ex2b, ex2c)
setGeneric("tre<-", function(x, value) standardGeneric("tre<-"))
#' @rdname assign-tre
#' @aliases tre<-,phyloseq,phylo-method
setMethod("tre<-", c("phyloseq", "phylo"), function(x, value){
	phyloseq(otuTable=x@otuTable, samData=x@samData, taxTab=x@taxTab, tre=value)
})
#' @rdname assign-tre
#' @aliases tre<-,phyloseq,phyloseq-method
setMethod("tre<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(otuTable=x@otuTable, samData=x@samData, taxTab=x@taxTab, tre=tre(value))
})
################################################################################
# # # # # assign.internal <- function(x, value, slot){
	# # # # # arglist <- splat.phyloseq.objects(x)
	# # # # # arglist[slot] <- value
	# # # # # do.call("phyloseq", arglist)
# # # # # }
# # # # # assign.internal(x, value, "tre")
# # # # # assign.internal(x, value, "otuTable")
# # # # # assign.internal(x, value, "")