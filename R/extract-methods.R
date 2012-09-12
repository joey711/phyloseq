################################################################################
# subsetting functions
# Without these, the default coerces to the base object (e.g. matrix or data.frame)
################################################################################
#' Extract parts of otu_table
#'
#' @export
#' @aliases [,otu_table-method
#' @rdname extract-methods
setMethod("[", "otu_table", function(x, i, j, ...){
	newx <- as(x, "matrix")[i, j, drop=FALSE]
	otu_table(newx, taxa_are_rows(x) )
})
################################################################################
#' extract parts of sample_data
#'
#' @export
#' @aliases [,sample_data-method
#' @rdname extract-methods
setMethod("[", "sample_data", function(x, i, j, ...){
	sample_data( data.frame(x)[i, j, drop=FALSE] )
})
################################################################################
#' extract parts of taxonomyTable
#'
#' @export
#' @aliases [,taxonomyTable-method
#' @rdname extract-methods
setMethod("[", "taxonomyTable", function(x, i, j, ...){
	tax_table( as(x, "matrix")[i, j, drop=FALSE] )
})
################################################################################
################################################################################
#' Generic extraction from higher-order object
#'
#' @export
#' @aliases [,phyloseq-method
#' @rdname extract-methods
setMethod("[", "phyloseq", function(x, i, j, ...){
	argslist <- list(...)
	# # # return(argslist)
	goodComponents <- c(get.component.classes(), names(get.component.classes()))
	goodComponents <- unique(goodComponents)
	goodComponents <- goodComponents[!goodComponents %in% c("phy_tree", "phylo")]

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
	
	# idiosyncracy, component class and slot name do not match for taxonomyTable/tax_table
	if( component == "taxonomyTable" ){
		component <- "tax_table"
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
