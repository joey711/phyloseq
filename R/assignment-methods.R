################################################################################
#' Assign a new OTU Table to \code{x}
#'
#' @usage otu_table(x) <- value
#'
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required).
#'  \code{\link{otu_table-class}}
#'  or
#'  \code{\link{phyloseq-class}}.
#'
#' @export
#' @docType methods
#' @rdname assign-otu_table
#' @aliases assign-otu_table otu_table<- otuTable<-
#'
#' @examples
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # OTU <- otu_table(GlobalPatterns)[1:100, ]
#' # otu_table(ex2b) <- OTU
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Relace otu_table by implying the component in context.
#' # ex2c <- GlobalPatterns
#' # otu_table(ex2c) <- ex2b
#' # identical(ex2a, ex2c)
setGeneric("otu_table<-", function(x, value) standardGeneric("otu_table<-"))
#' @rdname assign-otu_table
#' @aliases otu_table<-,phyloseq,otu_table-method
setMethod("otu_table<-", c("phyloseq", "otu_table"), function(x, value){
	phyloseq(value, x@sam_data, x@tax_table, x@phy_tree, x@refseq)
})
#' @rdname assign-otu_table
#' @aliases otu_table<-,otu_table,otu_table-method
setMethod("otu_table<-", c("otu_table", "otu_table"), function(x, value){ value })
#' @rdname assign-otu_table
#' @aliases otu_table<-,phyloseq,phyloseq-method
setMethod("otu_table<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(otu_table(value), x@sam_data, x@tax_table, x@phy_tree, x@refseq)
})
################################################################################
#' Manually change taxa_are_rows through assignment.
#'
#' The taxa_are_rows slot is a logical indicating the orientation of the
#' abundance table contained in object \code{x}.
#'
#' @usage taxa_are_rows(x) <- value
#'
#' @param x \code{\link{otu_table-class}} or \code{\link{phyloseq-class}}
#'
#' @param value A logical of length equal to 1. If \code{length(value) > 1}, 
#'  the additional elements will be ignored. Only the first element is assigned
#'  to the taxa_are_rows slot.
#'
#' @export
#' @docType methods
#' @rdname assign-taxa_are_rows
#' @aliases assign-taxa_are_rows taxa_are_rows<- speciesAreRows<-
#'
#' @examples #
#' # data(GlobalPatterns)
#' # taxa_are_rows(GlobalPatterns)
#' # taxa_are_rows(otu_table(GlobalPatterns))
setGeneric("taxa_are_rows<-", function(x, value){
	standardGeneric("taxa_are_rows<-")
})
#' @rdname assign-taxa_are_rows
#' @aliases taxa_are_rows<-,otu_table,logical-method
setMethod("taxa_are_rows<-", c("otu_table", "logical"), function(x, value){
	x@taxa_are_rows <- value[1]
	return(x)
})
#' @rdname assign-taxa_are_rows
#' @aliases taxa_are_rows<-,phyloseq,logical-method
setMethod("taxa_are_rows<-", c("phyloseq", "logical"), function(x, value){
	taxa_are_rows(otu_table(x)) <- value
	return(x)
})
################################################################################
#' Assign (new) sample_data to \code{x}
#'
#' This replaces the current \code{sample_data} component of \code{x} with 
#' \code{value}, if \code{value} is a \code{\link{sample_data-class}}. However,
#' if \code{value} is a \code{data.frame}, then \code{value} is first coerced to
#' a \code{\link{sample_data-class}}, and then assigned. Alternatively, if 
#' \code{value} is \code{\link{phyloseq-class}}, then the 
#' \code{\link{sample_data}} component will first be accessed from \code{value}
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
#' @usage sample_data(x) <- value
#'
#' @param x (Required). \code{\link{phyloseq-class}}. The object to modify.
#' @param value (Required). Either a \code{\link{sample_data-class}}, 
#'  a \code{data.frame} that can be coerced into \code{\link{sample_data-class}}, 
#'  or a \code{\link{phyloseq-class}} that contains a 
#'  suitable \code{sample_data} component to assign to \code{x}. If unsure,
#'  try \code{\link{sample_data}}\code{(value)}, which should return a 
#'  \code{\link{sample_data-class}} object without error.
#'
#' @return No return. This is an assignment statement.
#'
#' @export
#' @rdname assign-sample_data
#' @aliases assign-sample_data sample_data<- sam_data<- sampleData<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 10 samples in GlobalPatterns
#' # ex2a <- prune_samples(sample_names(GlobalPatterns)[1:10], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # SD <- sample_data(GlobalPatterns)[1:10, ]
#' # sample_data(ex2b) <- SD
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Example restoring the original sample_data component. ex2c lacks sample_data
#' # ex2c <- phyloseq(otu_table(GlobalPatterns), tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
#' # sample_data(ex2c) <- GlobalPatterns
#' # identical(ex2c, GlobalPatterns)
#' # # Can try on ex2b, but other components have only 10 samples. No change.
#' # sample_data(ex2b) <- GlobalPatterns
#' # identical(ex2a, ex2b) # still true.
"sample_data<-" <- function(x, value){
	if( !inherits(value, "sample_data") ){
		value <- sample_data(value)
	}
	phyloseq(x@otu_table, value, x@tax_table, x@phy_tree, x@refseq)
}
################################################################################
#' Assign a (new) Taxonomy Table to \code{x}
#'
#' @usage tax_table(x) <- value
#' 
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{taxonomyTable-class}}.
#'  Alternatively, \code{value} can be a \code{\link{phyloseq-class}} that has
#'  a \code{\link{tax_table}} component, or a \code{\link{matrix-class}}
#'  that can be coerced to a \code{\link{taxonomyTable-class}} with row indices
#'  that match at least some of the \code{\link{taxa_names}} of \code{x}.
#'
#' @export
#' @rdname assign-tax_table
#' @aliases assign-tax_table tax_table<- taxTab<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # TT <- tax_table(GlobalPatterns)[1:100, ]
#' # tax_table(ex2b) <- TT
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # 2 examples adding a tax_table component from phyloseq or matrix classes
#' # ex2c <- phyloseq(otu_table(ex2b), sample_data(ex2b), phy_tree(ex2b))
#' # tax_table(ex2c) <- ex2b
#' # identical(ex2a, ex2c)
#' # ex2c <- phyloseq(otu_table(ex2b), sample_data(ex2b), phy_tree(ex2b))
#' # tax_table(ex2c) <- as(tax_table(ex2b), "matrix")
#' # identical(ex2a, ex2c)
setGeneric("tax_table<-", function(x, value) standardGeneric("tax_table<-"))
#' @rdname assign-tax_table
#' @aliases tax_table<-,phyloseq,taxonomyTable-method
setMethod("tax_table<-", c("phyloseq", "taxonomyTable"), function(x, value){
	phyloseq(x@otu_table, x@sam_data, value, x@phy_tree, x@refseq)
})
#' @rdname assign-tax_table
#' @aliases tax_table<-,phyloseq,ANY-method
setMethod("tax_table<-", c("phyloseq", "ANY"), function(x, value){
	phyloseq(x@otu_table, x@sam_data, tax_table(value, FALSE), x@phy_tree, x@refseq)
})
#' @rdname assign-tax_table
#' @aliases tax_table<-,taxonomyTable,taxonomyTable-method
setMethod("tax_table<-", c("taxonomyTable", "taxonomyTable"), function(x, value){
	# Asign as-is.
	value
})
#' @rdname assign-tax_table
#' @aliases tax_table<-,taxonomyTable,ANY-method
setMethod("tax_table<-", c("taxonomyTable", "ANY"), function(x, value){
	tax_table(value, FALSE)
})
################################################################################
#' Assign a (new) phylogenetic tree to \code{x}
#'
#' @usage phy_tree(x) <- value
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{phylo-class}}, or \code{\link{phyloseq-class}}
#'
#' @export
#' @docType methods
#' @rdname assign-phy_tree
#' @aliases assign-phy_tree phy_tree<- tre<-
#' @examples #
#' # data(GlobalPatterns)
#' # # An example of pruning to just the first 100 taxa in GlobalPatterns.
#' # ex2a <- prune_species(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
#' # # The following 3 lines produces an ex2b that is equal to ex2a
#' # ex2b <- GlobalPatterns
#' # tree <- prune_species(taxa_names(GlobalPatterns)[1:100], phy_tree(GlobalPatterns))
#' # phy_tree(ex2b) <- tree
#' # identical(ex2a, ex2b)
#' # print(ex2b)
#' # # Example adding a phylo tree from phyloseq class
#' # ex2c <- phyloseq(otu_table(ex2b), sample_data(ex2b), tax_table(ex2b))
#' # phy_tree(ex2c) <- ex2b
#' # identical(ex2b, ex2c)
setGeneric("phy_tree<-", function(x, value) standardGeneric("phy_tree<-"))
#' @rdname assign-phy_tree
#' @aliases phy_tree<-,phyloseq,phylo-method
setMethod("phy_tree<-", c("phyloseq", "phylo"), function(x, value){
	phyloseq(x@otu_table, x@sam_data, x@tax_table, value, x@refseq)
})
#' @rdname assign-phy_tree
#' @aliases phy_tree<-,phyloseq,phyloseq-method
setMethod("phy_tree<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(x@otu_table, x@sam_data, x@tax_table, phy_tree(value), x@refseq)
})
################################################################################
#' Replace OTU identifier names
#'
#' @usage taxa_names(x) <- value
#'
#' @param x (Required). An object defined by the \code{\link{phyloseq-package}}
#' 	that describes OTUs in some way.
#' @param value (Required). A character vector 
#'  to replace the current \code{\link{taxa_names}}.
#'
#' @export
#' @docType methods
#' @rdname assign-taxa_names
#' @aliases assign-taxa_names taxa_names<-
#'
#' @examples
#' data("esophagus")
#' taxa_names(esophagus)
#' # plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
#' taxa_names(esophagus) <- paste("OTU-", taxa_names(esophagus), sep="")
#' taxa_names(esophagus)
#' # plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
#' ## non-characters are first coerced to characters.
#' taxa_names(esophagus) <- 1:ntaxa(esophagus)
#' taxa_names(esophagus)
#' # plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
#' ## Cannot assign non-unique or differently-lengthed name vectors. Error.
#' # taxa_names(esophagus) <- sample(c(TRUE, FALSE), ntaxa(esophagus), TRUE)
#' # taxa_names(esophagus) <- sample(taxa_names(esophagus), ntaxa(esophagus)-5, FALSE)
setGeneric("taxa_names<-", function(x, value){
	if( anyDuplicated(value) ){
		stop("taxa_names<-: You are attempting to assign duplicated taxa_names")
	}
	standardGeneric("taxa_names<-")
})
# Attempt to coerce value to a character vector. Remaining methods will require it.
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,ANY,ANY-method
setMethod("taxa_names<-", c("ANY", "ANY"), function(x, value){
  taxa_names(x) <- as(value, "character")
  return(x)
})
# value is now character, but no specific method for first argumet
# return x unchanged.
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,ANY,character-method
setMethod("taxa_names<-", c("ANY", "character"), function(x, value){
  return(x)
})
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,otu_table,character-method
setMethod("taxa_names<-", c("otu_table", "character"), function(x, value){
  if( taxa_are_rows(x) ){
    rownames(x) <- value
  } else {
    colnames(x) <- value
  }
  return(x)
})
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,taxonomyTable,character-method
setMethod("taxa_names<-", c("taxonomyTable", "character"), function(x, value){
  rownames(x) <- value
  return(x)
})
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,phylo,character-method
setMethod("taxa_names<-", c("phylo", "character"), function(x, value){
  x$tip.label <- value
  return(x)
})
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,XStringSet,character-method
setMethod("taxa_names<-", c("XStringSet", "character"), function(x, value){
  names(x) <- value
  return(x)
})
#' @rdname assign-taxa_names
#' @aliases taxa_names<-,phyloseq,character-method
setMethod("taxa_names<-", c("phyloseq", "character"), function(x, value){
  # dispatch on components
  taxa_names(x@otu_table) <- value
  taxa_names(x@phy_tree)  <- value
  taxa_names(x@tax_table) <- value
  taxa_names(x@refseq)    <- value
  return(x)
})
################################################################################
################################################################################
#' Replace OTU identifier names
#'
#' @usage sample_names(x) <- value
#'
#' @param x (Required). An object defined by the \code{\link{phyloseq-package}}
#' 	that describes OTUs in some way.
#' @param value (Required). A character vector 
#'  to replace the current \code{\link{sample_names}}.
#'
#' @export
#' @docType methods
#' @rdname assign-sample_names
#' @aliases assign-sample_names sample_names<-
#'
#' @examples
#' data("esophagus")
#' sample_names(esophagus)
#' # plot_tree(esophagus, color="sample_names", ladderize="left")
#' sample_names(esophagus) <- paste("Sa-", sample_names(esophagus), sep="")
#' sample_names(esophagus)
#' # plot_tree(esophagus, color="sample_names", ladderize="left") 
#' ## non-characters are first coerced to characters.
#' sample_names(esophagus) <- 1:nsamples(esophagus)
#' sample_names(esophagus)
#' # plot_tree(esophagus, color="sample_names", ladderize="left") 
#' ## Cannot assign non-unique or differently-lengthed name vectors. Error.
#' # sample_names(esophagus) <- sample(c(TRUE, FALSE), nsamples(esophagus), TRUE)
#' # sample_names(esophagus) <- sample(sample_names(esophagus), nsamples(esophagus)-1, FALSE)
setGeneric("sample_names<-", function(x, value){
	if( anyDuplicated(value) ){
		stop("sample_names<-: You are attempting to assign duplicated sample_names")
	}
	standardGeneric("sample_names<-")
})
# Attempt to coerce value to a character vector. Remaining methods will require it.
#' @rdname assign-sample_names
#' @aliases sample_names<-,ANY,ANY-method
setMethod("sample_names<-", c("ANY", "ANY"), function(x, value){
	sample_names(x) <- as(value, "character")
	return(x)
})
# value is now character, but no specific method for first argumet
# return x unchanged.
#' @rdname assign-sample_names
#' @aliases sample_names<-,ANY,character-method
setMethod("sample_names<-", c("ANY", "character"), function(x, value){
	return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,otu_table,character-method
setMethod("sample_names<-", c("otu_table", "character"), function(x, value){
	if( taxa_are_rows(x) ){
		colnames(x) <- value
	} else {
		rownames(x) <- value
	}
	return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,sample_data,character-method
setMethod("sample_names<-", c("sample_data", "character"), function(x, value){
	rownames(x) <- value
	return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,phyloseq,character-method
setMethod("sample_names<-", c("phyloseq", "character"), function(x, value){
	# dispatch on components
	sample_names(x@otu_table) <- value
	sample_names(x@sam_data)  <- value
	return(x)
})
################################################################################