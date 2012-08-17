################################################################################
#' Assign a new OTU Table to \code{x}
#'
#' @usage otu_table(x) <- value
#'
#' @param x (Required). \code{\link{phyloseq-class}}
#' @param value (Required). \code{\link{otu_table-class}} or \code{\link{phyloseq-class}}.
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
	phyloseq(otu_table=value, sam_data=x@sam_data, tax_table=x@tax_table, phy_tree=x@phy_tree)
})
#' @rdname assign-otu_table
#' @aliases otu_table<-,otu_table,otu_table-method
setMethod("otu_table<-", c("otu_table", "otu_table"), function(x, value){ value })
#' @rdname assign-otu_table
#' @aliases otu_table<-,phyloseq,phyloseq-method
setMethod("otu_table<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(otu_table=otu_table(value), sam_data=x@sam_data, tax_table=x@tax_table, phy_tree=x@phy_tree)
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
	if(class(value) != "sample_data"){value <- sample_data(value)}
	phyloseq(otu_table=x@otu_table, sam_data=value, tax_table=x@tax_table, phy_tree=x@phy_tree)
}
#' @export
#' @rdname assign-sample_data
#' @aliases assign-sample_data sample_data<- sam_data<-
#' @usage sam_data(x) <- value
"sam_data<-" <- function(x, value){
	if(class(value) != "sample_data"){value <- sample_data(value)}
	phyloseq(otu_table=x@otu_table, sam_data=value, tax_table=x@tax_table, phy_tree=x@phy_tree)
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
	phyloseq(otu_table=x@otu_table, sam_data=x@sam_data, tax_table=value, phy_tree=x@phy_tree)
})
#' @rdname assign-tax_table
#' @aliases tax_table<-,phyloseq,ANY-method
setMethod("tax_table<-", c("phyloseq", "ANY"), function(x, value){
	phyloseq(otu_table=x@otu_table, sam_data=x@sam_data, tax_table=tax_table(value, FALSE), phy_tree=x@phy_tree)
})
#' @rdname assign-tax_table
#' @aliases tax_table<-,taxonomyTable,taxonomyTable-method
setMethod("tax_table<-", c("taxonomyTable", "taxonomyTable"), function(x, value){
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
	phyloseq(otu_table=x@otu_table, sam_data=x@sam_data, tax_table=x@tax_table, phy_tree=value)
})
#' @rdname assign-phy_tree
#' @aliases phy_tree<-,phyloseq,phyloseq-method
setMethod("phy_tree<-", c("phyloseq", "phyloseq"), function(x, value){
	phyloseq(otu_table=x@otu_table, sam_data=x@sam_data, tax_table=x@tax_table, phy_tree=phy_tree(value))
})
################################################################################
