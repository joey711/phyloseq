################################################################################
#' Assign to otuTable an object/value.
#'
#' @usage otuTable(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @docType methods
#' @rdname assign-otuTable
#' @aliases assign-otuTable otuTable<-
#'
#' @examples #
setGeneric("otuTable<-", function(x, value) standardGeneric("otuTable<-"))
#' @rdname assign-otuTable
#' @aliases otuTable<-,phyloseq,otuTable-method
setMethod("otuTable<-", c("phyloseq", "otuTable"), function(x, value){
	phyloseq(otuTable=value, sampleMap=x@sampleMap, taxTab=x@taxTab, tre=x@tre)
})
################################################################################
#' Manually change speciesAreRows through assignment.
#'
#' The speciesAreRows slot is a logical indicating the orientation of the
#' abundance table contained in object \code{x}.
#'
#' @usage speciesarerows(x) <- value
#'
#' @param x An otuTable-class or higher-order object that contains an otuTable.
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
#' Assign to sampleMap an object/value.
#'
#' @usage sampleMap(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-sampleMap
#' @aliases assign-sampleMap sampleMap<-
#' @examples #
"sampleMap<-" <- function(x, value){
	phyloseq(otuTable=x@otuTable, sampleMap=value, taxTab=x@taxTab, tre=x@tre)
}
################################################################################
#' Assign to taxTab an object/value.
#'
#' @usage taxTab(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-taxTab
#' @aliases assign-taxTab taxTab<-
#' @examples #
"taxTab<-" <- function(x, value){
	phyloseq(otuTable=x@otuTable, sampleMap=x@sampleMap, taxTab=value, tre=x@tre)
}
################################################################################
#' Assign to tre an object/value.
#'
#' This will automatically convert "phylo"-class trees (ape package) to
#' "phylo4"-class trees (phylobase package).
#'
#' @usage tre(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @docType methods
#' @rdname assign-tre
#' @aliases assign-tre tre<-
#' @examples #
setGeneric("tre<-", function(x, value) standardGeneric("tre<-"))
#' @rdname assign-tre
#' @aliases tre<-,phyloseq,phylo-method
setMethod("tre<-", c("phyloseq", "phylo"), function(x, value){
	phyloseq(otuTable=x@otuTable, sampleMap=x@sampleMap, taxTab=x@taxTab, tre=as(value, "phylo4"))
})
#' @rdname assign-tre
#' @aliases tre<-,phyloseq,phylo4-method
setMethod("tre<-", c("phyloseq", "phylo4"), function(x, value){
	phyloseq(otuTable=x@otuTable, sampleMap=x@sampleMap, taxTab=x@taxTab, tre=value)
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