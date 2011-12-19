################################################################################
#' The S4 class for storing taxa-abundance information.
#'
#' Because orientation of these tables can vary by method, the orientation is
#' defined explicitly in the \code{speciesAreRows} slot (a logical).
#' The \code{otuTable} class inherits the \code{\link{matrix}} class to store
#' abundance values.
#' Various standard subset and assignment nomenclature has been extended to apply
#' to the \code{otuTable} class, including square-bracket, \code{\link{t}}, etc.
#'
#' \describe{
#'    \item{speciesAreRows}{
#'		A single logical specifying the orientation of the abundance table.
#'    }
#'
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#'  }
#' @name otuTable-class
#' @rdname otuTable-class
#' @exportClass otuTable
setClass("otuTable", representation(speciesAreRows="logical"), contains = "matrix")
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
#' @name sampleData-class
#' @rdname sampleData-class
#' @exportClass sampleData
setClass("sampleData", contains="data.frame")
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
setClassUnion("otuTableOrNULL", c("otuTable", "NULL"))
#' @keywords internal
setClassUnion("sampleDataOrNULL", c("sampleData", "NULL"))
#' @keywords internal
setClassUnion("taxonomyTableOrNULL", c("taxonomyTable", "NULL"))
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
################################################################################
# The actual phyloseq master class with all 4 slots. This is akin to 
# the otuSamTaxTree class of previous versions, but 
# with the possibility of empty (NULL) slots and an explicit prototype 
# for slots to be NULL if they are not provided at instantiation.
################################################################################
#' The main experiment-level class for phyloseq data
#'
#' Contains all component classes: 
#' \code{\link{otuTable-class}},
#' \code{\link{sampleData-class}},
#' \code{\link{taxonomyTable-class}} (\code{"taxTab"} slot), and
#' \code{\link[ape]{phylo}}-class (\code{"tre"} slot). There are several advantages
#' to storing your phylogenetic sequencing experiment as an instance of the
#' phyloseq class, not the least of which is that it is easy to return to the
#' data later and feel confident that the different data types ``belong'' to
#' one another. Furthermore, the \code{\link{phyloseq}} constructor ensures that
#' the different data components have compatible indices (e.g. species and samples),
#' and performs the necessary trimming automatically when you create your
#' ``experiment-level'' object. Downstream analyses are aware of which data
#' classes they require -- and where to find them -- often making your 
#' \code{phyloseq-class} object the only data argument to analysis and plotting
#' functions (although there are many options and parameter arguments waiting
#' for you). 
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
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{samData}{ a single object of class sampleData.}
#'    \item{taxTab}{ a single object of class taxonomyTable.}
#'    \item{tre}{ a single object of class phylo, from the package ape}
#' }
#' @seealso The constructor, \code{\link{phyloseq}}, 
#'  the merger \code{\link{merge_phyloseq}}, and also the component 
#'  constructor/accessors \code{\link{otuTable}}, \code{\link{sampleData}},
#'  \code{\link{taxTab}}, and \code{\link{tre}}.
#' 
#' @name phyloseq-class
#' @rdname phyloseq-class
#' @exportClass phyloseq
setClass(Class="phyloseq", 
	representation=representation(
		otuTable="otuTableOrNULL",
		taxTab="taxonomyTableOrNULL",
		samData="sampleDataOrNULL",
		tre="phyloOrNULL"),
	prototype=prototype(otuTable=NULL, taxTab=NULL, samData=NULL, tre=NULL)
)
################################################################################