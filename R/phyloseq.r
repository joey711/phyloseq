################################################################################
#' Build objects for the phyloseq package.
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
	
	splatlist <- sapply(arglist, splat.phyloseq.objects)

	# Name hash:
	slot_hash <- c("otu", "Sam", "Tax", "Tree")
	names(slot_hash) <- c("otuTable", "sampleMap", "taxTab", "tre")
	newClassName <- slot_hash[names(splatlist)]
	classNameOrder <- sapply(newClassName, 
		function(i, slot_hash){which(slot_hash==i)},
	slot_hash)
	names(classNameOrder) <- newClassName
	classNameOrder <- sort(classNameOrder)
	newClassName   <- paste(names(classNameOrder), sep="", collapse="")

	## Need to determine which new() type to call.
	# First some quality-control checks:
	if( length(splatlist) > 4){
		stop("phyloseq()-ERROR: Too many components provided\n")
	} else if( length(names(splatlist)) > length(unique(names(splatlist))) ){
		stop("phyloseq()-ERROR: Only one of each component type allowed.\n",
		"For merging multiple tables of the same class, try merge_phyloseq(...)\n")
	} else if( length(splatlist) == 1){
		return(arglist[[1]])
	} else if( !"otuTable" %in% names(splatlist) ){
		stop("phyloseq()-ERROR: Argument list must include an otuTable.\n")
	# And finally, instantiate and return the H.O. object.
	} else {
		return( do.call("new", c(list(Class=newClassName), splatlist) ) )
	}
}
################################################################################