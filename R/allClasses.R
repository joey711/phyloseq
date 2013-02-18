################################################################################
#' The S4 class for storing taxa-abundance information.
#'
#' Because orientation of these tables can vary by method, the orientation is
#' defined explicitly in the \code{taxa_are_rows} slot (a logical).
#' The \code{otu_table} class inherits the \code{\link{matrix}} class to store
#' abundance values.
#' Various standard subset and assignment nomenclature has been extended to apply
#' to the \code{otu_table} class, including square-bracket, \code{\link{t}}, etc.
#'
#' \describe{
#'    \item{taxa_are_rows}{
#'		A single logical specifying the orientation of the abundance table.
#'    }
#'
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#'  }
#' @name otu_table-class
#' @rdname otu_table-class
#' @exportClass otu_table
setClass("otu_table", representation(taxa_are_rows="logical"), contains = "matrix")
################################################################################
#' The S4 for storing sample variables.
#'
#' Row indices represent samples, while column indices represent experimental
#' categories, variables (and so forth) that describe the samples.
#'
#' \describe{
#'
#'    \item{.Data}{data-frame data, inherited from the data.frame class.}
#' 
#'    \item{row.names}{
#'	     Also inherited from the data.frame class;
#'       it should contain the sample names.
#'    }
#' 
#'    \item{names}{Inherited from the data.frame class.}
#' 
#'  }
#' 
#' @name sample_data-class
#' @rdname sample_data-class
#' @exportClass sample_data
setClass("sample_data", contains="data.frame")
################################################################################
#' An S4 class that holds taxonomic classification data as a character
#' matrix.
#'
#' Row indices represent taxa, columns represent taxonomic classifiers.
#' 
#' \describe{
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#' }
#'
#' @name taxonomyTable-class
#' @rdname taxonomyTable-class
#' @exportClass taxonomyTable
setClass("taxonomyTable", contains = "matrix")
################################################################################
#' An S4 copy of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree. It is used throught ape.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
#' @name phylo-class
#' @rdname phylo-class
#' @exportClass phylo
setOldClass("phylo")
################################################################################
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.
################################################################################
#' @keywords internal
setClassUnion("otu_tableOrNULL", c("otu_table", "NULL"))
#' @keywords internal
setClassUnion("sample_dataOrNULL", c("sample_data", "NULL"))
#' @keywords internal
setClassUnion("taxonomyTableOrNULL", c("taxonomyTable", "NULL"))
#' @import ape
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
#' @import Biostrings
#' @keywords internal
setClassUnion("XStringSetOrNULL", c("XStringSet", "NULL"))
################################################################################
#' The main experiment-level class for phyloseq data
#'
#' Contains all currently-supported component data classes: 
#' \code{\link{otu_table-class}},
#' \code{\link{sample_data-class}},
#' \code{\link{taxonomyTable-class}} (\code{"tax_table"} slot),
#' \code{\link[ape]{phylo}}-class (\code{"phy_tree"} slot),
#' and the \code{\link[Biostrings]{XStringSet-class}} (\code{"refseq"} slot).
#' There are several advantages
#' to storing your phylogenetic sequencing experiment as an instance of the
#' phyloseq class, not the least of which is that it is easy to return to the
#' data later and feel confident that the different data types ``belong'' to
#' one another. Furthermore, the \code{\link{phyloseq}} constructor ensures that
#' the different data components have compatible indices (e.g. OTUs and samples),
#' and performs the necessary trimming automatically when you create your
#' ``experiment-level'' object. Downstream analyses are aware of which data
#' classes they require -- and where to find them -- often making your 
#' \code{phyloseq-class} object the only data argument required for analysis and plotting
#' functions (although there are many options and parameter arguments available
#' to you). 
#'
#' In the case of missing component data, the slots are set to \code{NULL}. As
#' soon as a \code{phyloseq-class} object is to be updated with new component
#' data (previously missing/\code{NULL} or not), the indices of all components
#' are re-checked for compatibility and trimmed if necessary. This is to ensure
#' by design that components describe the same taxa/samples, and also that these
#' trimming/validity checks do not need to be repeated in downstream analyses.
#' 
#' slots:
#' \describe{
#'    \item{otu_table}{a single object of class otu_table.}
#'    \item{sam_data}{ a single object of class sample_data.}
#'    \item{tax_table}{ a single object of class taxonomyTable.}
#'    \item{phy_tree}{ a single object of the \code{\link[ape]{phylo}}-class, from the ape package.}
#'    \item{refseq}{ a biological sequence set object of a class that
#'         inherits from the \code{\link[Biostrings]{XStringSet-class}}, from the Biostrings package.}
#' }
#' @seealso
#'  The constructor, \code{\link{phyloseq}}, 
#'  the merger \code{\link{merge_phyloseq}}, and also the component 
#'  constructor/accessors \code{\link{otu_table}}, \code{\link{sample_data}},
#'  \code{\link{tax_table}}, \code{\link{phy_tree}}, and \code{\link{refseq}}.
#' 
#' @import ape
#' @import Biostrings
#' @name phyloseq-class
#' @rdname phyloseq-class
#' @exportClass phyloseq
setClass(Class="phyloseq", 
	representation=representation(
		otu_table="otu_tableOrNULL",
		tax_table="taxonomyTableOrNULL",
		sam_data="sample_dataOrNULL",
		phy_tree="phyloOrNULL",
		refseq = "XStringSetOrNULL"),
	prototype=prototype(otu_table=NULL, tax_table=NULL, sam_data=NULL, phy_tree=NULL, refseq=NULL)
)
################################################################################
