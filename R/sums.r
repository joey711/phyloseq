################################################################################
#' Returns the total number of individuals observed from each species.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables
#'
#' @usage speciesSums(x)
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
speciesSums <- function(x){
	x <- otuTable(x)
	if( speciesAreRows(x) ){
		rowSums(x)
	} else {
		colSums(x)
	}
}
speciessums <- speciesSums
################################################################################
#' Returns the total number of individuals observed from each sample.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables.
#'
#' @usage sampleSums(x)
#'
#' @param x Any phyloseq-package object that is or contains an otuTable.
#' 
#' @return the total number of individuals present in each sample. A named 
#'  integer vector with length equal to the number of samples
#'  in the table, name indicated the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @seealso speciesSums rowSums colSums sum
#' @export
#' @examples #
sampleSums <- function(x){
	x <- otuTable(x)
	if( speciesAreRows(x) ){
		colSums(x)
	} else {
		rowSums(x)
	}
}
samplesums <- sampleSums
################################################################################
