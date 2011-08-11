################################################################################
#' Returns the total number of individuals observed from each species.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables
#'
#' @param x Any phyloseq-package object that is or contains an otuTable.
#' 
#' @return A named integer vector with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the sum of
#'  all individuals observed for each taxa in \code{x}.
#'
#' @seealso sampleSums rowSums colSums
#' @export
#' @examples #
setGeneric("speciesSums", function(x) standardGeneric("speciesSums"))
setMethod("speciesSums", "otuTable", function(x){
	if( speciesAreRows(x) ){
		species_sums = rowSums(x)
	} else {
		species_sums = colSums(x)
	}
	return(species_sums)
}) 
setMethod("speciesSums", "otuTree",  function(x){ speciesSums(otuTable(x)) })
setMethod("speciesSums", "phyloseq", function(x){ speciesSums(otuTable(x)) })
speciessums <- speciesSums
################################################################################
################################################################################
