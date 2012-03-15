############################################################################
#' @aliases show,otuTable-method
#' @rdname show-methods
setMethod("show", "otuTable", function(object){
	# print otuTable (always there).
	cat(paste("OTU Table:          [", nspecies(object), " species and ", 
        nsamples(object), " samples]", sep = ""), fill = TRUE)	
    if( speciesAreRows(object) ){
    	cat("                     species are rows", fill=TRUE)
    } else { 
    	cat("                     species are columns", fill=TRUE)
    }
	show(as(object, "matrix"))
})
############################################################################
#' @aliases show,sampleData-method
#' @rdname show-methods
setMethod("show", "sampleData", function(object){
	cat(paste("Sample Data:         [", dim(sampleData(object))[1], " samples by ", 
		dim(sampleData(object))[2], 
		" sample variables]:", sep = ""),
	fill = TRUE)
	show(as(object, "data.frame"))
})
############################################################################
#' @aliases show,taxonomyTable-method
#' @rdname show-methods
setMethod("show", "taxonomyTable", function(object){
    cat(paste("Taxonomy Table:     [", dim(object)[1], " species by ", 
		dim(object)[2], 
		" taxonomic ranks]:", sep = ""),
	fill = TRUE)	
	show(as(object, "matrix"))
})
############################################################################
#' method extensions to show for phyloseq objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @export
#' @aliases show,phyloseq-method
#' @docType methods
#' @rdname show-methods
#' @examples
#' # data(GlobalPatterns)
#' # show(GlobalPatterns)
#' # GlobalPatterns
setMethod("show", "phyloseq", function(object){
	cat("phyloseq-class experiment-level object", fill=TRUE)
	# print otuTable (always there).
	cat(paste("OTU Table:          [", nspecies(otuTable(object)), " species and ", 
        nsamples(otuTable(object)), " samples]", sep = ""), fill = TRUE)	
    if( speciesAreRows(otuTable(object)) ){
    	cat("                     species are rows", fill=TRUE)
    } else { 
    	cat("                     species are columns", fill=TRUE)
    } 

	# print Sample Data if there
	if(!is.null(sampleData(object, FALSE))){
        cat(paste("Sample Map:         [", dim(sampleData(object))[1], " samples by ", 
	        dim(sampleData(object))[2], 
            " sample variables]:", sep = ""), fill = TRUE)
	}

	# print tax Tab if there	
	if(!is.null(taxTab(object, FALSE))){
        cat(paste("Taxonomy Table:     [", dim(taxTab(object))[1], " species by ", 
	        dim(taxTab(object))[2], 
            " taxonomic ranks]:", sep = ""), fill = TRUE)
	}
	
	# print tree if there
	if(!is.null(tre(object, FALSE))){
        cat(paste("Phylogenetic Tree:  [", nspecies(tre(object)), " tips and ", 
	        tre(object)$Nnode,
            " internal nodes]", sep = ""),
        	fill = TRUE
        )
		if( is.rooted(tre(object)) ){
			cat("                     rooted", fill=TRUE)
		} else {
			cat("                     unrooted", fill=TRUE)
		}        
	}
})
############################################################################