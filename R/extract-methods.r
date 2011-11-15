################################################################################
# subsetting functions
# Without these, the default coerces to the base object (e.g. matrix or data.frame)
################################################################################
#' Extract parts of otuTable
#'
#' @export
#' @aliases [,otuTable,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "otuTable", function(x,i,j,...){
	newx <- callNextMethod(x@.Data,i,j,drop=FALSE,...)
	new("otuTable", x@speciesAreRows, newx)
})
################################################################################
#' extract parts of sampleMap
#'
#' @export
#' @aliases [,sampleMap,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "sampleMap", function(x,i,j,...){
	new("sampleMap", callNextMethod(data.frame(x),i,j,drop=FALSE,...))
})
################################################################################
#' extract parts of taxonomyTable
#'
#' @export
#' @aliases [,taxonomyTable,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "taxonomyTable", function(x,i,j,...){
	new("taxonomyTable", callNextMethod(x@.Data,i,j,drop=FALSE,...))
})
################################################################################
################################################################################
#' Generic extraction from higher-order object
#'
#' @export
#' @aliases [,phyloseqFather,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "phyloseqFather", function(x, i, j, ...){
	argslist <- list(...)
	# # # return(argslist)
	goodComponents <- c(get.component.classes(), names(get.component.classes()))
	goodComponents <- unique(goodComponents)
	goodComponents <- goodComponents[!goodComponents %in%
						c("tre", "old-tre", "phylo", "phylo4")]

	# If there is an argument labeled "component", use that. Else, search among unlabled
	if(any( names(argslist)=="component" )){
		component <- argslist[["component"]]
		if( !component %in% goodComponents ){
			stop("component not valid. Please select from\n",
				"  the non-tree components returned by get.component.classes()")
		}
	# Else, if there are unlabeled arguments, check those
	} else if( length(argslist) > 0  & 
				any(sapply(argslist, class)=="character")  & 
				any(argslist %in% goodComponents) ){
		# If there is a component class supplied in arg list, take the first
		component <- argslist[[ which(argslist %in% goodComponents)[1] ]]
	# else, punt, return warning()
	} else {
		stop("Attempt to extract from complex S4-phyloseq object without\n", 
		"specifying the component from which you would like to extract.\n",
		"Try get.component.classes() to see possible values to provide.")
	}
	
	# idiosyncracy, component class and slot name do not match for taxonomyTable/taxTab
	if( component == "taxonomyTable" ){
		component <- "taxTab"
	}
	newx    <- access(x, component)
	
	# Nested-if to protect against missing arguments. Can only miss one at a time.
	if( missing(i) ){
		return( newx[, j] )
	} else if( missing(j) ){
		return( newx[i, ] )
	} else {
		return( newx[i, j] )
	}
})
################################################################################
################################################################################