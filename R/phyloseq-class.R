################################################################################
#' Build phyloseq-class objects from their components.
#'
#' \code{phyloseq()} is a constructor method, This is the main method
#' suggested for constructing an experiment-level (\code{\link{phyloseq-class}})
#' object from its component data 
#' (component data classes: \code{\link{otu_table-class}}, \code{\link{sample_data-class}}, 
#'  \code{\link{taxonomyTable-class}}, \code{\link{phylo-class}}).
#'
#' @usage phyloseq(...)
#'
#' @param ... One or more component objects among the set of classes
#'  defined by the phyloseq package, as well as \code{phylo}-class
#'  (defined by the \code{\link{ape-package}}). Each argument should be a different class.
#'  For combining multiple components of the same class, or multiple phyloseq-class
#'  objects, use the \code{\link{merge_phyloseq}} function. Unlike in earlier
#'  versions, the arguments to phyloseq do not need to be named, and the order
#'  of the arguments does not matter.
#'
#' @return The class of the returned object depends on the argument 
#'  class(es). For an experiment-level object, two or more component data objects
#'  must be provided.
#'  Otherwise, if a single component-class
#'  is provided, it is simply returned as-is. 
#'  The order of arguments does not matter. 
#'
#' @seealso \code{\link{merge_phyloseq}}
#' @export
#' @examples #
#' # # data(GlobalPatterns)
#' # # GP <- GlobalPatterns
#' # # phyloseq(sample_data(GP), otu_table(GP))
#' # # phyloseq(otu_table(GP), phy_tree(GP))
#' # # phyloseq(tax_table(GP), otu_table(GP))
#' # # phyloseq(phy_tree(GP), otu_table(GP), sample_data(GP))
#' # # phyloseq(otu_table(GP), tax_table(GP), sample_data(GP))
#' # # phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))
phyloseq <- function(...){
	
	arglist <- list(...)
	
	# Remove names from arglist. Will replace them based on their class
	names(arglist) <- NULL

	# ignore all but component data classes.
	arglist  <- arglist[sapply(arglist, class) %in% phyloseq:::get.component.classes()]
	
	# Make the name-replaced, splatted list
	splatlist <- sapply(arglist, phyloseq:::splat.phyloseq.objects)

	## Need to determine which new() type to call.
	# First some quality-control checks:
	if( length(splatlist) > 4){
		stop("Too many components provided\n")
	} else if( length(names(splatlist)) > length(unique(names(splatlist))) ){
		stop("Only one of each component type allowed.\n",
		"For merging multiple objects of the same class, try merge_phyloseq(...)\n")
	} else if( length(splatlist) == 1){
		return(arglist[[1]])
	# Instantiate the phyloseq-class object, ps.
	} else {
		ps <- do.call("new", c(list(Class="phyloseq"), splatlist) )
	}
	# Verify there is more than one component that describes species before attempting to reconcile.
	if( sum(!sapply(lapply(phyloseq:::splat.phyloseq.objects(ps), species.names), is.null)) >= 2 ){	
		ps <- phyloseq:::reconcile_species(ps)
	}
	# Verify there is more than one component that describes samples before attempting to reconcile.
	if( sum(!sapply(lapply(phyloseq:::splat.phyloseq.objects(ps), sample.names), is.null)) >= 2 ){
		ps <- phyloseq:::reconcile_samples(ps)		
	}
	# ENFORCE CONSISTENT ORDER OF TAXA INDICES.
	# If there is a phylogenetic tree included, re-order the otu_table based 
	# according to the order of taxa-names on the tree, and optionally for
	# the taxonomyTable, if present.
	if( !is.null(phy_tree(ps, FALSE)) ){
		otu <- as(otu_table(ps), "matrix")
		# Re-order the matrix order matches tree
		if( taxa_are_rows(ps) ){
			otu <- otu[species.names(phy_tree(ps)), ]
		} else {
			otu <- otu[, species.names(phy_tree(ps))]
		}
		ps@otu_table <- otu_table(otu, taxa_are_rows(ps))
		
		# If there is a taxonomyTable, re-order that too.
		if( !is.null(tax_table(ps, FALSE)) ){
			tax <- as(tax_table(ps), "matrix")
			tax <- tax[species.names(phy_tree(ps)), ]
			ps@tax_table <- tax_table(tax)
		}
	}

	# ENFORCE CONSISTENT ORDER OF SAMPLE INDICES
	# Other errors have been detected for when sample indices do not match.
	# check first that ps has sample_data
	if( !is.null(sample_data(ps, FALSE)) ){
		if( !all(sample.names(otu_table(ps)) == rownames(sample_data(ps))) ){
			# Reorder the sample_data rows so that they match the otu_table order.
			ps@sam_data <- sample_data(ps)[sample.names(otu_table(ps)), ]
		}		
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
#' element is named by the corresponding slot name in the phyloseq-class.
#'
#' @keywords internal
#'
#' @examples #
#' #get.component.classes()
get.component.classes <- function(){
	# define classes vector
	component.classes <- c("otu_table", "sample_data", "phylo", "taxonomyTable")
	# the names of component.classes needs to be the slot names to match getSlots / splat
	names(component.classes) <- c("otu_table", "sam_data", "phy_tree", "tax_table")	
	return(component.classes)
}
################################################################################
#' Convert \code{\link{phyloseq-class}} into a named list of its non-empty components.
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
	# Check if class of x is among the component classes (not phyloseq-class)
	if( class(x) %in% component.classes ){
		splatx <- list(x)
		names(splatx) <- names(component.classes)[component.classes==class(x)]
	} else if( class(x) == "phyloseq" ){ 
		slotnames <- names(getSlots("phyloseq"))
		allslots  <- sapply(slotnames, function(i, x){access(x, i, FALSE)}, x)
		splatx    <- allslots[!sapply(allslots, is.null)]
	} else {
		return(NULL)
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
#' the \code{\link{phyloseq-class}}.
#' 
#' @seealso merge_phyloseq
#' @export
#' @examples #
#'  data(GlobalPatterns)
#'  getslots.phyloseq(GlobalPatterns)
#'  data(esophagus)
#'  getslots.phyloseq(esophagus)
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
#' data component. Unlike the component-specific accessors (e.g. \code{\link{otu_table}},
#' or \code{\link{phy_tree}}),
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
#' ## data(GlobalPatterns)
#' ## access(GlobalPatterns, "tax_table")
#' ## access(GlobalPatterns, "phy_tree")
#' ## access(otu_table(GlobalPatterns), "otu_table")
#' ## # Should return NULL:
#' ## access(otu_table(GlobalPatterns), "sample_data")
#' ## access(otuTree(GlobalPatterns), "sample_data")
#' ## access(otuSam(GlobalPatterns), "phy_tree")
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
#' @param x (Required). A \code{\link{phyloseq-class}} object
#'  that contains 2 or more components
#'  that in-turn describe species/taxa.
#'
#' @return Returns a character vector of only those species that are present in
#'  all species-describing components of \code{x}.
#'
#' @seealso \code{\link{reconcile_species}}, \code{\link{Reduce}}
#' @keywords internal
#' @examples #
#' ## data(GlobalPatterns)
#' ## head(intersect_species(GlobalPatterns), 10)
intersect_species <- function(x){
	component_list  <- splat.phyloseq.objects(x)
	doesnt_have_species <- which( getslots.phyloseq(x) %in% c("sam_data") )
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
#'  \code{\link{otu_table}}, \code{\link{tax_table}}, \code{\link{phy_tree}}.
#'
#' @return A trimmed version of the argument, \code{x}, in which each component
#'  describes exactly the same set of species/taxa. Class of return should match
#'  argument, \code{x}.
#'
#' @seealso \code{\link{reconcile_samples}}, \code{\link{Reduce}}
#' @rdname reconcile_species-methods
#' @keywords internal
#' @examples #
#' ## data(GlobalPatterns)
#' ## head(intersect_species(GlobalPatterns), 10)
#' ## reconcile_species(GlobalPatterns)
#' # # data(phylocom)
#' # # tree <- phylocom$phylo
#' # # OTU  <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
#' # # ex3  <- phyloseq(OTU, tree)
#' # # reconcile_species(ex3)
#' # # intersect_species(ex3)
setGeneric("reconcile_species", function(x) standardGeneric("reconcile_species"))
################################################################################
#' @aliases reconcile_species,phyloseq-method
#' @rdname reconcile_species-methods
setMethod("reconcile_species", signature("phyloseq"), function(x){
	# Make species the intersection of all species in the components
	species <- intersect_species(x)
	# prune_species(species, x) already checks if species and species.names(x)
	# sets are equal, no need to re-check.
	x       <- prune_species(species, x)
	return(x)
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
#' ## data(GlobalPatterns)
#' ## reconcile_samples(GlobalPatterns)
reconcile_samples <- function(x){
	# prevent infinite recursion issues by checking if intersection already satisfied
	if( setequal(sample.names(sample_data(x)), sample.names(otu_table(x))) ){
		return(x)
	} else {
		samples <- intersect(sample.names(sample_data(x)), sample.names(otu_table(x)))		
		x@sam_data <- prune_samples(samples, x@sam_data)
		x@otu_table  <- prune_samples(samples, x@otu_table)
		return(x)
	}
}
################################################################################
