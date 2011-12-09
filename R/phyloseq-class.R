################################################################################
#' Build phyloseq-class objects from their components.
#'
#' \code{phyloseq()} is a constructor method, This is the main method
#' suggested for constructing \code{phyloseq} higher-order
#' objects from their components.
#'
#' @usage phyloseq(...)
#'
#' @param ... One or more component objects among the set of classes
#'  defined by the phyloseq package, as well as \code{phylo4} objects
#'  (defined by the phylobase package). Each argument should be a different class.
#'  For combining multiple components of the same class, or multiple higher-order
#'  classes, use the \code{\link{merge_phyloseq}} function. Unlike in earlier
#'  versions, the arguments to phyloseq do not need to be named, and the order
#'  of the arguments does not matter.
#'
#' @return The class of the returned object depends on the argument 
#'  class(es). To construct a H.O. object, two or more component data objects
#'  must be provided in the argument list.
#'  Otherwise, the order of arguments does not matter. If a single component-class
#'  object is provided, it is simply returned as-is. 
#'
#' @seealso \code{\link{merge_phyloseq}}
#' @export
#' @examples #
#' # # data(ex1)
#' # # phyloseq(sampleMap(ex1), otuTable(ex1))
#' # # phyloseq(otuTable(ex1), tre(ex1))
#' # # phyloseq(taxTab(ex1), otuTable(ex1))
#' # # phyloseq(tre(ex1), otuTable(ex1), sampleMap(ex1))
#' # # phyloseq(otuTable(ex1), taxTab(ex1), sampleMap(ex1))
#' # # phyloseq(otuTable(ex1), tre(ex1), taxTab(ex1), sampleMap(ex1))
phyloseq <- function(...){
	
	arglist <- list(...)

	## Check for "phylo" class objects. If present, convert them to "phylo4"
	if( any(sapply(arglist, class)=="phylo") ){
		for( i in which(sapply(arglist, class)=="phylo") ){
			arglist[[i]] <- as(arglist[[i]], "phylo4")
		}
	}
	
	# Remove names from arglist. Will replace them based on their class
	names(arglist) <- NULL
	
	# Make the name-replaced, splatted list, vetted by splat.phlyoseq.objects
	splatlist <- sapply(arglist, splat.phyloseq.objects)

	## Need to determine which new() type to call.
	# First some quality-control checks:
	if( length(splatlist) > 4){
		stop("phyloseq()-ERROR: Too many components provided\n")
	} else if( length(names(splatlist)) > length(unique(names(splatlist))) ){
		stop("phyloseq()-ERROR: Only one of each component type allowed.\n",
		"For merging multiple tables of the same class, try merge_phyloseq(...)\n")
	} else if( length(splatlist) == 1){
		return(arglist[[1]])
	} else if( !("otuTable" %in% sapply(splatlist, class)) ){
		stop("phyloseq()-ERROR: Argument list must include an otuTable.\n")
	# Instantiate the phyloseq-class object, ps.
	} else {
		ps <- do.call("new", c(list(Class="phyloseq"), splatlist) )
	}
	# Verify there is more than one component that describes species before attempting to reconcile.
	if( sum(!sapply(lapply(splat.phyloseq.objects(ps), species.names), is.null)) >= 2 ){	
		ps <- reconcile_species(ps)
	}
	# Verify there is more than one component that describes samples before attempting to reconcile.
	if( sum(!sapply(lapply(splat.phyloseq.objects(ps), sample.names), is.null)) >= 2 ){
		ps <- reconcile_samples(ps)		
	}
	return(ps)
}
################################################################################
################################################################################
#' Show the component objects classes and slot names.
#'
#' There are no arguments to this function. It returns a named character
#' when called, which can then be used for tests of component data types, etc.
#'
#' @usage get.component.classes()
#' 
#' @return a character vector of the component objects classes, where each 
#' element is named by the corresponding slot name in the higher-order 
#' phyloseq objects (objects containing more than 1 phyloseq data object).
#' @keywords internal
#' @examples #
#' #get.component.classes()
get.component.classes <- function(){
	# define classes vector
	component.classes <- c("otuTable", "sampleMap", "phylo4", "phylo", "taxonomyTable")
	# the names of component.classes needs to be the slot names to match getSlots / splat
	names(component.classes) <- c("otuTable", "sampleMap", "tre", "old-tre", "taxTab")	
	return(component.classes)
}
################################################################################
#' Convert phyloseq-class into a named list of its non-empty components.
#'
#' This is used in internal handling functions, and one of its key features
#' is that the names in the returned-list match the slot-names, which is useful
#' for constructing calls with language-computing functions like \code{\link{do.call}}.
#' Another useful aspect is that it only returns the contents of non-empty slots.
#' In general, this should only be used by phyloseq-package developers. Standard
#' users should not need or use this function, and should use the accessors and 
#' other tools that leave the multi-component object in one piece.
#'
#' @usage splat.phyloseq.objects(x)
#'
#' @param x A \code{\link{phyloseq-class}} object. Alternatively, a component
#'  data object will work, resulting in named list of length 1. 
#' 
#' @return A named list, where each element is a component object that was contained 
#' in the argument, \code{x}. Each element is named according to its slot-name in
#' the phyloseq-object from which it is derived. 
#' If \code{x} is already a component data object,
#' then a list of length (1) is returned, also named.
#'
#' @seealso merge_phyloseq
#' @keywords internal
#' @examples #
splat.phyloseq.objects <- function(x){
	component.classes <- get.component.classes()
	# Check if class of x is among the component classes (not H.O.)
	if( class(x) %in% component.classes ){
		splatx <- list(x)
		names(splatx) <- names(component.classes)[component.classes==class(x)]
	} else { 
		slotnames <- names(getSlots(class(x)))
		allslots  <- sapply(slotnames, function(i, x){access(x, i, FALSE)}, x)
		splatx    <- allslots[!sapply(allslots, is.null)]
	}
	return(splatx)
}
################################################################################
#' Return the non-empty slot names of a phyloseq object.
#'
#' Like \code{\link{getSlots}}, but returns the class name if argument 
#' is component data object.
#' 
#' @usage getslots.phyloseq(physeq)
#'
#' @param physeq A \code{\link{phyloseq-class}} object. If \code{physeq} is a component
#'  data class, then just returns the class of \code{physeq}.
#' 
#' @return identical to getSlots. A named character vector of the slot classes
#' of a particular S4 class, where each element is named by the slot name it
#' represents. If \code{physeq} is a component data object,
#' then a vector of length (1) is returned, named according to its slot name in
#' the higher-order objects.
#' 
#' @seealso merge_phyloseq
#' @export
#' @examples #
getslots.phyloseq <- function(physeq){
	# Check if class of physeq is among the component classes (not phyloseq-class)
	component.classes <- get.component.classes()	
	if( class(physeq) %in% component.classes ){
		slotsx        <- as.character(class(physeq))
		names(slotsx) <- names(component.classes)[component.classes==class(physeq)]
	} else {
		# Make sure to return only the names of non-empty slots of physeq
		splatx <- splat.phyloseq.objects(physeq)
		slotsx <- names(splatx)
	}
	return(slotsx)
}
################################################################################
#' Universal slot accessor function for phyloseq-class.
#'
#' This function is used internally by many accessors and in 
#' many functions/methods that need to access a particular type of component data.
#' If something is wrong, or the slot is missing, the expected behavior is that
#' this function will return NULL. Thus, the output can be tested by 
#' \code{\link{is.null}} as verification of the presence of a particular 
#' data component. Unlike the component-specific accessors (e.g. \code{\link{otuTable}},
#' or \code{\link{tre}}),
#' the default behavior is not to stop with an error if the desired slot is empty.
#' In all cases this is controlled by the \code{errorIfNULL} argument, which can
#' be set to \code{TRUE} if an error is desired. 
#'
#' @usage access(physeq, slot, errorIfNULL=FALSE)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#'
#' @param slot (Required). A character string indicating the slot (not data class)
#'  of the component data type that is desired.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{FALSE}. 
#'
#' @return Returns the component object specified by the argument \code{slot}. 
#'  Returns NULL if slot does not exist. Returns \code{physeq} as-is 
#'  if it is a component class that already matches the slot name.
#'
#' @seealso \code{\link{getslots.phyloseq}}, \code{\link{merge_phyloseq}}
#' @export
#' @examples #
#' ## data(ex1)
#' ## access(ex1, "taxTab")
#' ## access(ex1, "tre")
#' ## access(otuTable(ex1), "otuTable")
#' ## # Should return NULL:
#' ## access(otuTable(ex1), "sampleMap")
#' ## access(otuTree(ex1), "sampleMap")
#' ## access(otuSam(ex1), "tre")
access <- function(physeq, slot, errorIfNULL=FALSE){
	component.classes <- get.component.classes()
	# Check if class of x is among the component classes (not H.O.)
	if( class(physeq) %in% component.classes ){
		# if slot-name matches physeq, return physeq as-is.
		if( component.classes[slot] == class(physeq) ){
			out <- physeq
		} else {
			out <- NULL
		}
	} else if(!slot %in% slotNames(physeq) ){
		out <- NULL
	} else {
		out <- eval(parse(text=paste("physeq@", slot, sep=""))) 
	}
	# Test if you should error upon the emptiness of the slot being accessed
	if( errorIfNULL & is.null(out) ){
		stop(slot, " slot is empty.")
	}
	return(out)
}
################################################################################
################################################################################
#' Returns the intersection of species for the components of x
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a phyloseq-object have exactly the same taxa/species.
#' It relies heavily on the \code{\link{Reduce}} function to determine the 
#' strictly common species.
#'
#' @usage intersect_species(x)
#'
#' @param x An object of the phyloseq package that contains 2 or more components
#'  data tables that in-turn describe species/taxa. E.g. An otuTree object, or
#'  an otuTax object.
#'
#' @return Returns a character vector of only those species that are present in
#'  all species-describing components of \code{x}.
#'
#' @seealso \code{\link{reconcile_species}}, \code{\link{Reduce}}
#' @keywords internal
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
#' @usage reconcile_species(x)
#'
#' @param x (Required). A phyloseq object. Only meaningful if \code{x} has at
#'  least two non-empty slots of the following slots that describe species:
#'  \code{\link{otuTable}}, \code{\link{taxTab}}, \code{\link{tre}}.
#'
#' @return A trimmed version of the argument, \code{x}, in which each component
#'  describes exactly the same set of species/taxa. Class of return should match
#'  argument, \code{x}.
#'
#' @seealso \code{\link{reconcile_samples}}, \code{\link{Reduce}}
#' @rdname reconcile_species-methods
#' @keywords internal
#' @examples #
#' ## data(ex1)
#' ## head(intersect_species(ex1), 10)
#' ## reconcile_species(ex1)
#' # # data(phylocom)
#' # # tree <- phylocom$phylo
#' # # OTU  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # # ex3  <- phyloseq(OTU, tree)
#' # # reconcile_species(ex3)
#' # # intersect_species(ex3)
setGeneric("reconcile_species", function(x) standardGeneric("reconcile_species"))
################################################################################
#' @aliases reconcile_species,phyloseq-method
#' @rdname reconcile_species-methods
setMethod("reconcile_species", signature("phyloseq"), function(x){
	# # # # # reconcile_species <- function(x){
	# Make species the intersection of all species in the components
	species    <- intersect_species(x)
	# Make allspecies the union of all species in the components,
	# the default for species.names() if argument is a phyloseq-class
	allspecies <- species.names(x)
	# prevent infinite recursion issues and unecessary mucking around
	# by checking if intersection and union are already the same (same species, all)
	if( setequal(species, allspecies) ){
		return(x)
	} else {
		x <- prune_species(species, x)
		return(x)
	}
})
################################################################################
#' Keep only sample-indices common to all components.
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a phyloseq-object describe exactly the same samples.
#' In expected cases, a user will not need to invoke this
#' function, because phyloseq objects are reconciled during instantiation by
#' default.
#'
#' @usage reconcile_samples(x)
#'
#' @param x An instance of phyloseq-class that contains 2 or more component
#'  data tables that in-turn describe samples. 
#'
#' @return A trimmed version of the argument, \code{x}, in which each component
#'  describes exactly the same set of samples. Class of \code{x} should be
#'  unchanged.
#'
#' @seealso \code{\link{reconcile_species}}
#' @keywords internal
#' @examples #
#' ## data(ex1)
#' ## reconcile_samples(ex1)
reconcile_samples <- function(x){
	# prevent infinite recursion issues by checking if intersection already satisfied
	if( setequal(sample.names(sampleMap(x)), sample.names(otuTable(x))) ){
		return(x)
	} else {
		samples <- intersect(sample.names(sampleMap(x)), sample.names(otuTable(x)))		
		x@sampleMap <- prune_samples(samples, x@sampleMap)
		x@otuTable  <- prune_samples(samples, x@otuTable)
		return(x)
	}
}
################################################################################
