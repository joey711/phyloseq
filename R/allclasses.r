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
#' Keystone S4 virtual class inherited by all other complex phyloseq classes.
#'
#' All other complex classes inherit by virtue of the fact that they 
#' contain an otuTable as one of their data components.
#' 
#' Slots:
#' \describe{
#'    \item{otuTable}{ Contains an otu abundance table object. See otuTable class.}
#'  }
#' @name phyloseqFather-class
#' @rdname phyloseqFather-class
#' @exportClass phyloseqFather
setClass(Class="phyloseqFather", 
	representation=representation(otuTable="otuTable", "VIRTUAL")
)
################################################################################
#' Define the otuSam class.
#'
#' Contains otuTable and sampleMap slots.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{a single object of class sampleMap.}
#' }
#' @name otuSam-class
#' @rdname otuSam-class
#' @exportClass otuSam
setClass(Class="otuSam", 
	representation=representation(sampleMap="sampleMap"),
	contains="phyloseqFather"
) 
# OLD
#setClass("phyloseq", representation(otuTable="otuTable", sampleMap="sampleMap"))
################################################################################
#' Define the otuTree class.
#'
#' Contains otuTable, phylo (tre) slots.
#' 
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{tre}{a single object of class phylo, from the package ape.}
#'   }
#'
#' @name otuTree-class
#' @rdname otuTree-class
#' @exportClass otuTree
setClass(Class="otuTree", representation=representation(tre="phylo4"), contains="phyloseqFather")
# OLD
#setClass("otuTree",  representation(otuTable="otuTable", tre="phylo"))
################################################################################
#' Define the otuTax class.
#'
#' Contains otuTable, phylo (tre) slots.
#' 
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{taxTab}{a single object of class taxonomyTable.}
#'   }
#'
#' @name otuTax-class
#' @rdname otuTax-class
#' @exportClass otuTax
setClass(Class="otuTax", representation=representation(taxTab="taxonomyTable"), contains="phyloseqFather")
# OLD
#setClass(Class="otuTax", representation=representation(otuTable="otuTable", taxTab="taxonomyTable"))
################################################################################
#' Define the otuSamTree class.
#'
#' Contains otuTable, sampleMap, phylo (tre) slots. Inherits the otuSam and
#' otuTree classes.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable}
#'    \item{sampleMap}{a single object of class sampleMap.}
#'    \item{tre}{a single object of class phylo, from the package ape}
#' }
#' @name otuSamTree-class
#' @rdname otuSamTree-class
#' @exportClass otuSamTree
setClass("otuSamTree", contains=c("otuSam", "otuTree"))
################################################################################
#' Define the otuSamTax class.
#'
#' Inherits the otuSam class.
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{a single object of class sampleMap.}
#'    \item{taxTab}{a single object of class taxonomyTable.}
#' }
#' @name otuSamTax-class
#' @rdname otuSamTax-class
#' @exportClass otuSamTax
setClass("otuSamTax",  contains=c("otuSam", "otuTax"))
# OLD
#setClass("otuSamTax", representation(taxTab="taxonomyTable"), contains="otuSam")	
################################################################################
#' Define the otuSamTaxTree class.
#'
#' Contains all (current) component classes: otuTable, sampleMap,
#' taxonomyTable (taxTab), and phylo (tre). 
#' Inherits otuSamTax and otuSamTree.
#' 
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{sampleMap}{ a single object of class sampleMap.}
#'    \item{taxTab}{ a single object of class taxonomyTable.}
#'    \item{tre}{ a single object of class phylo, from the package ape}
#' }
#' @name otuSamTaxTree-class
#' @rdname otuSamTaxTree-class
#' @exportClass otuSamTaxTree
setClass("otuSamTaxTree", contains=c("otuSamTax", "otuSamTree"))	
