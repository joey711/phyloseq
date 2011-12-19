################################################################################
#' Assign a new OTU Table to \code{x}
#'
#' @usage otuTable(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{otuTable-class}}
#'
#' @export
#' @docType methods
#' @rdname assign-otuTable
#' @aliases assign-otuTable otuTable<-
#'
#' @examples
#' # data(ex1)
#' # # An example of pruning to just the first 100 taxa in ex1.
#' # ex2a <- prune_species(species.names(ex1)[1:100], ex1)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- ex1
#' # OTU <- otuTable(ex1)[1:100, ]
#' # otuTable(ex2b) <- OTU
#' # identical(ex2a, ex2b)
#' # print(ex2b)
setGeneric("otuTable<-", function(x, value) standardGeneric("otuTable<-"))
#' @rdname assign-otuTable
#' @aliases otuTable<-,phyloseq,otuTable-method
setMethod("otuTable<-", c("phyloseq", "otuTable"), function(x, value){
	phyloseq(otuTable=value, samData=x@samData, taxTab=x@taxTab, tre=x@tre)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuTable,otuTable-method
setMethod("otuTable<-", c("otuTable", "otuTable"), function(x, value){ value })
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
#' # data(ex1)
#' # speciesarerows(ex1)
#' # speciesarerows(otuTable(ex1))
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
#' Assign a new sampleData to \code{x}
#'
#' @usage sampleData(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{sampleData-class}}
#'
#' @export
#' @rdname assign-sampleData
#' @aliases assign-sampleData sampleData<- samData<-
#' @examples #
#' # data(ex1)
#' # # An example of pruning to just the first 10 samples in ex1
#' # ex2a <- prune_samples(sample.names(ex1)[1:10], ex1)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- ex1
#' # SD <- sampleData(ex1)[1:10, ]
#' # sampleData(ex2b) <- SD
#' # identical(ex2a, ex2b)
#' # print(ex2b)
"sampleData<-" <- function(x, value){
	phyloseq(otuTable=x@otuTable, samData=value, taxTab=x@taxTab, tre=x@tre)
}
#' @export
#' @rdname assign-sampleData
#' @aliases assign-sampleData sampleData<- samData<-
#' @usage samData(x) <- value
"samData<-" <- function(x, value){
	phyloseq(otuTable=x@otuTable, samData=value, taxTab=x@taxTab, tre=x@tre)
}
################################################################################
#' Assign a new Taxonomy Table to \code{x}
#'
#' @usage taxTab(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{taxonomyTable-class}}
#'
#' @export
#' @rdname assign-taxTab
#' @aliases assign-taxTab taxTab<-
#' @examples #
#' # data(ex1)
#' # # An example of pruning to just the first 100 taxa in ex1.
#' # ex2a <- prune_species(species.names(ex1)[1:100], ex1)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- ex1
#' # TT <- taxTab(ex1)[1:100, ]
#' # taxTab(ex2b) <- TT
#' # identical(ex2a, ex2b)
#' # print(ex2b)
"taxTab<-" <- function(x, value){
	phyloseq(otuTable=x@otuTable, samData=x@samData, taxTab=value, tre=x@tre)
}
################################################################################
#' Assign to tre an object/value.
#'
#' @usage tre(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{phylo-class}}
#'
#' @export
#' @docType methods
#' @rdname assign-tre
#' @aliases assign-tre tre<-
#' @examples #
#' # data(ex1)
#' # # An example of pruning to just the first 100 taxa in ex1.
#' # ex2a <- prune_species(species.names(ex1)[1:100], ex1)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- ex1
#' # tree <- prune_species(species.names(ex1)[1:100], tre(ex1))
#' # tre(ex2b) <- tree
#' # identical(ex2a, ex2b)
#' # print(ex2b)
setGeneric("tre<-", function(x, value) standardGeneric("tre<-"))
#' @rdname assign-tre
#' @aliases tre<-,phyloseq,phylo-method
setMethod("tre<-", c("phyloseq", "phylo"), function(x, value){
	phyloseq(otuTable=x@otuTable, samData=x@samData, taxTab=x@taxTab, tre=value)
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