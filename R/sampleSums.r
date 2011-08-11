################################################################################
#' Returns the total number of individuals observed from each sample.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables.
#'
#' @param x Any phyloseq-package object that is or contains an otuTable.
#' 
#' @return the total number of individuals present in each sample. A named 
#'  integer vector with length equal to the number of samples
#'  in the table, name indicated the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @seealso speciesSums rowSums colSums sum
#' @keywords sum sums samples
#' @export
#' @examples #
setGeneric("sampleSums", function(x) standardGeneric("sampleSums"))
setMethod("sampleSums", "otuTable", function(x){
	if( speciesAreRows(x) ){
		sample_sums = colSums(x)
	} else {
		sample_sums = rowSums(x)
	}
	return(sample_sums)
}) 
setMethod("sampleSums", "otuTree",  function(x){ sampleSums(otuTable(x)) })
setMethod("sampleSums", "phyloseq", function(x){ sampleSums(otuTable(x)) })
################################################################################
################################################################################
