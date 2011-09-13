################################################################################
#' Returns the intersection of species for the components of x
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a phyloseq-object have exactly the same taxa/species.
#' It relies heavily on the \code{\link{Reduce}} function to determine the 
#' strictly common species.
#'
#' @param x An object of the phyloseq package that contains 2 or more components
#'  data tables that in-turn describe species/taxa. E.g. An otuTree object, or
#'  an otuTax object.
#'
#' @return Returns a character vector of only those species that are present in
#'  all species-describing components of \code{x}.
#'
#' @seealso \code{\link{reconcile_species}}, \code{\link{Reduce}}
#' @export
#' @examples #
#' ## data(ex1)
#' ## head(intersect_species(ex1), 10)
intersect_species <- function(x){
	component_list  <- splat.phyloseq.objects(x)
	doesnt_have_species <- which( getslots.phyloseq(x) %in% c("sampleMap") )
	if( length(doesnt_have_species) > 0 ){
		species_vectors <- lapply(component_list[-doesnt_have_species], species.names)		
	} else {
		species_vectors <- lapply(component_list, species.names)		
	}
	return( Reduce("intersect", species_vectors) )
}
################################################################################
#' Keep only species-indices common to all components.
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a phyloseq-object have exactly the same taxa/species.
#' It relies heavily on the \code{\link{prune_species}} S4 methods to perform the
#' actual trimming. In expected cases, a user will not need to invoke this
#' function, because phyloseq objects are reconciled during instantiation by
#' default.
#'
#' @param x An object of the phyloseq package that contains 2 or more components
#'  data tables that in-turn describe species/taxa. E.g. An otuTree object, or
#'  an otuTax object.
#'
#' @return A trimmed version of the argument, \code{x}, in which each component
#'  describes exactly the same set of species/taxa. Class of \code{x} should be
#'  unchanged.
#'
#' @seealso \code{\link{reconcile_samples}}, \code{\link{Reduce}}
#' @export
#' @examples #
#' ## data(ex1)
#' ## head(intersect_species(ex1), 10)
#' ## reconcile_species(ex1)
reconcile_species <- function(x){
	species <- intersect_species(x)
	# prevent infinite recursion issues by checking if intersection already satisfied
	if( setequal(species, species.names(x)) ){
		return(x)
	} else {
		slots2reconcile <- names(getslots.phyloseq(x))
		for( i in slots2reconcile ){
			eval(parse(text=paste("x@", i, "<- prune_species(species, x@", i, ")", sep="")))
		}
		return(x)
	}
}
################################################################################
#' Keep only sample-indices common to all components.
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a phyloseq-object describe exactly the same samples.
#' In expected cases, a user will not need to invoke this
#' function, because phyloseq objects are reconciled during instantiation by
#' default.
#'
#' @param x An object of the phyloseq package that contains 2 or more components
#'  data tables that in-turn describe samples. At present this is limited to 
#'  the otuSam class (otuTable and SampleMap) and its children.
#'
#' @return A trimmed version of the argument, \code{x}, in which each component
#'  describes exactly the same set of samples. Class of \code{x} should be
#'  unchanged.
#'
#' @seealso \code{\link{reconcile_species}}
#' @export
#' @examples #
#' ## data(ex1)
#' ## reconcile_samples(ex1)
reconcile_samples <- function(x){
	samples <- intersect(rownames(x@sampleMap), sample.names(x@otuTable))
	# prevent infinite recursion issues by checking if intersection already satisfied
	if( setequal(samples, sample.names(x)) ){
		return(x)
	} else {
		x@sampleMap <- prune_samples(samples, x@sampleMap)
		x@otuTable  <- prune_samples(samples, x@otuTable)
		return(x)
	}
}
################################################################################
