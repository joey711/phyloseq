################################################################################
#' Prune unwanted species / taxa from a phylogenetic object.
#' 
#' An S4 Generic method for removing (pruning) unwanted spceies from phylogenetic
#' objects, including phylo and phylo4 trees, as well as native phyloseq package
#' objects. This is particularly useful for pruning objects of the 
#' \code{\link{otuTable}} and the higher order objects that contain it as a component.
#' Note, the \code{phylo} class version is adapted from \code{picante::prune.samples}
#'
#' @param species A character vector of the species in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept species are TRUE, and length
#' is equal to the number of species in object x. If \code{species} is a named
#' logical, the species retained is based on those names. Make sure they are
#' compatible with the \code{species.names} of the object you are modifying (\code{x}). 
#'
#' @param x A phylogenetic object, including \code{phylo} and \code{phylo4} trees,
#' as well as native phyloseq package objects that represent taxa / species
#' (that is, otuTable, taxonomyTable, phyloseq, otuTree, and their superclasses).
#'
#' @return The class of the object returned by \code{prune_species} matches
#' the class of the argument, \code{x}.
#'
#' @keywords prune trim subset species OTU taxa
#' @export
#' @examples #
#' ## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' ## f1  <- filterfunSample(topk(2))
#' ## wh1 <- genefilterSample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_species(wh1, testOTU)
#' ## prune_species(wh2, testOTU)
#' ## 
#' ## taxtab1 <- taxTab(matrix("abc", 5, 5))
#' ## prune_species(wh1, taxtab1)
#' ## prune_species(wh2, taxtab1)
setGeneric("prune_species", function(species, x) standardGeneric("prune_species"))
setMethod("prune_species", signature("character", "phylo"), function(species, x){
	require("picante")
	trimTaxa <- setdiff(x$tip.label, species)
	if( length(trimTaxa) > 0 ){
		drop.tip(x, trimTaxa)
	} else x
})
setMethod("prune_species", signature("character", "otuTree"), function(species, x){
	tre(x) <- prune_species( species, tre(x) )
	return(x)
})

setMethod("prune_species", signature("character", "phylo4"), function(species, x){
	trimTaxa <- setdiff(tipLabels(x), species)
	if( length(trimTaxa) > 0 ){
		subset(x, tips.exclude=trimTaxa) # drop.tip(x, trimTaxa)
	} else x
})
setMethod("prune_species", signature("character", "otuTree4"), function(species, x){
	tre(x) <- prune_species( species, tre(x) )
	return(x)
})

# otuTable related methods
setMethod("prune_species", signature("logical", "otuTable"), function(species, x){
	if( is.null(names(species)) ){
		species.names(x) <- species.names(x)[species]
	} else {
		species.names(x) <- names(species)[species]
	}
	return(x)
})
setMethod("prune_species", signature("logical", "phyloseq"), function(species, x){
	otuTable(x) <- prune_species(species, otuTable(x))
	return(x)
})
setMethod("prune_species", signature("logical", "otuTree"), function(species, x){
	otuTable(x) <- prune_species(species, otuTable(x))
	return(x)
})

# taxonomyTable 
setMethod("prune_species", signature("logical", "taxonomyTable"), function(species, x){
	if( is.null(names(species)) ){
		x <- x[species, ]
	} else {
		x <- x[names(species)[species], ]
	}
	return(x)
})
################################################################################
