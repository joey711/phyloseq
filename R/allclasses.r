################################################################################
#' The S4 class that holds taxa-abundance information.
#'
#' Because orientation of these tables can vary by method, the orientation is
#' defined explicitly in the \code{speciesAreRows} slot (a logical).
#' The \code{otuTable} class inherits the \code{\link{matrix}} class to store
#' abundance values.
#' Various standard subset and assignment nomenclature has been extended to apply
#' to the \code{otuTable} class, including square-bracket, \code{\link{t}}, etc.
#'
#' \describe{
#'    \item{speciesAreRows}{A single logical specifying the orientation of the 
#' abundance table.}
#'
#'    \item{nspecies}{A single positive integer specifying the number of distinct 
#' species/taxa.}
#' 
#'    \item{nsamples}{ A single positive integer specifying the number of distinct
#' samples. Both \code{nsamples} and \code{nspecies} should not be assigned directly by user, but
#' are assessed internally during initialization of \code{otuTable} objects.}
#' 
#'    \item{species.names}{ A character vector of the names of the taxa / species,
#' so that they can be accessed, modified, subsetted, without the need for
#' the user to consider the orientation of the abundance table.}
#' 
#'    \item{sample.names}{A character vector of the names of the samples, so that they
#' can be accessed, modified, subsetted, without the need for the user to consider
#' the orientation of the abundance table. Both \code{sample.names} and
#' \code{species.names} can
#' be assigned, e.g. \code{species.names(object)<-}, however, this has the effect
#' of internally replacing the otuTable to be a subset that matches the new name 
#' vector, rather than replacing the index names in-place.}
#' 
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#'  }
#' @name otuTable-class
#' @rdname otuTable-class
#' @exportClass otuTable
setClass("otuTable",
	representation(speciesAreRows="logical",
		nspecies="integer",
		nsamples="integer",
		species.names="character",
		sample.names="character"),
	contains = "matrix"
)
################################################################################
#' The S4 object in the phyloseq package that holds sample data as a data.frame.
#'
#' Row indices represent samples, while column indices represent experimental
#' categories, variables, etc. that describe the samples.
#'
#' \describe{
#'    \item{nsamples}{ A single positive integer specifying the number of distinct samples.
#' \code{nsamples} should not be assigned directly by user, but
#' is assessed internally during initialization of \code{sampleMap} objects.}
#'
#'    \item{.Data}{data-frame data, inherited from the data.frame class.}
#' 
#'    \item{row.names}{Also inherited from the data.frame class;
#'  it should contain the sample names.}
#' 
#'    \item{names}{Inherited from the data.frame class.}
#' 
#'    \item{nsamples}{A single positive integer, indicating the number of samples
#'     (rows) represented by an object. Consistent with otuTable objects. 
#'      should not be assigned directly by user, but
#'      is assessed internally during instantiation of \code{sampleMap} objects.}
#'  }
#' 
#' @name sampleMap-class
#' @rdname sampleMap-class
#' @exportClass sampleMap
setClass("sampleMap", representation(nsamples="integer"), contains="data.frame")
################################################################################
############################################################################
#' An S4 class that holds taxonomic classification data as a character
#' matrix.
#'
#' Row indices represent taxa, columns represent taxonomic classifiers.
#' 
#' \describe{
#'    \item{nspecies}{A single positive integer specifying the number of 
#' distinct species/taxa. \code{nspecies} should not be assigned directly by user, but
#' is assessed internally during instantiation of \code{taxonomyTable}
#' objects.}
#'
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#' }
#'
#' @name taxonomyTable-class
#' @rdname taxonomyTable-class
#' @exportClass taxonomyTable
setClass("taxonomyTable", representation(nspecies="integer"), contains = "matrix")
############################################################################
############################################################################
############################################################################
# Define the higher-order classes in the phyloseq package.
# 
# 
############################################################################
#' Define the phyloseq class.
#'
#' Contains otuTable and sampleMap slots.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{a single object of class sampleMap.}
#' }
#' @name phyloseq-class
#' @rdname phyloseq-class
#' @exportClass phyloseq
setClass("phyloseq", representation(otuTable="otuTable", sampleMap="sampleMap"))
########################################
#' Define the otuTree class.
#'
#' Contains otuTable, phylo (tre) slots.
#' 
#' \describe{
#'    \item{otuTable}{a single object of class otuTable}
#'    \item{tre}{a single object of class phylo, from the package ape}
#'   }
#'
#' @name otuTree-class
#' @rdname otuTree-class
#' @exportClass otuTree
setClass("otuTree",  representation(otuTable="otuTable", tre="phylo"))				
########################################
#' Define the phyloseqTree class.
#'
#' Contains otuTable, sampleMap, phylo (tre) slots. Inherits the phyloseq and
#' otuTree classes.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable}
#'    \item{sampleMap}{a single object of class sampleMap.}
#'    \item{tre}{a single object of class phylo, from the package ape}
#' }
#' @name phyloseqTree-class
#' @rdname phyloseqTree-class
#' @exportClass phyloseqTree
setClass("phyloseqTree", contains=c("phyloseq", "otuTree"))
########################################
#' Define the phyloseqTax class.
#'
#' Inherits the phyloseq class.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{a single object of class sampleMap.}
#'    \item{taxTab}{a single object of class taxonomyTable.}
#' }
#' @name phyloseqTax-class
#' @rdname phyloseqTax-class
#' @exportClass phyloseqTax
setClass("phyloseqTax", representation(taxTab="taxonomyTable"), contains="phyloseq")	
########################################
#' Define the phyloseqTaxTree class.
#'
#' Contains all (current) component classes: otuTable, sampleMap,
#' taxonomyTable (taxTab), and phylo (tre). 
#' Inherits phyloseqTax and phyloseqTree.
#' 
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{ a single object of class sampleMap.}
#'    \item{taxTab}{ a single object of class taxonomyTable.}
#'    \item{tre}{ a single object of class phylo, from the package ape}
#' }
#' @name phyloseqTaxTree-class
#' @rdname phyloseqTaxTree-class
#' @exportClass phyloseqTaxTree
setClass("phyloseqTaxTree", contains=c("phyloseqTax", "phyloseqTree"))	
########################################
#' Define the otuTree4 class.
#'
#' Identical to otuTree class, but uses \code{phylo4} class for tree slot
#' instead of \code{phylo}.
#' \describe{
#'    \item{otuTable}{ a single object of class otuTable}
#'    \item{tre}{ a single object of class phylo4, imported from phylobase.}
#' }
#' @exportClass otuTree4
#'
#' @name otuTree4-class
#' @rdname otuTree4-class
setClass("otuTree4",  representation(otuTable="otuTable", tre="phylo4"))
### # # @importClassesFrom phylobase phylo4

