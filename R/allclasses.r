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
#'    \item{speciesAreRows}{
#'		A single logical specifying the orientation of the abundance table.
#'    }
#'
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#'  }
#' @name otuTable-class
#' @rdname otuTable-class
#' @exportClass otuTable
setClass("otuTable", representation(speciesAreRows="logical", contains = "matrix") )
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
setClass("sampleMap", contains="data.frame")
################################################################################
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
setClass("taxonomyTable", contains = "matrix")
################################################################################

################################################################################