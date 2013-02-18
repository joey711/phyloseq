############################################################################
#' @aliases show,otu_table-method
#' @rdname show-methods
setMethod("show", "otu_table", function(object){
	# print otu_table (always there).
	cat(paste("OTU Table:          [", ntaxa(object), " taxa and ", 
        nsamples(object), " samples]", sep = ""), fill = TRUE)	
    if( taxa_are_rows(object) ){
    	cat("                     taxa are rows", fill=TRUE)
    } else { 
    	cat("                     taxa are columns", fill=TRUE)
    }
	show(as(object, "matrix"))
})
############################################################################
#' @aliases show,sample_data-method
#' @rdname show-methods
setMethod("show", "sample_data", function(object){
	cat(paste("Sample Data:        [", dim(sample_data(object))[1], " samples by ", 
		dim(sample_data(object))[2], 
		" sample variables]:", sep = ""),
	fill = TRUE)
	show(as(object, "data.frame"))
})
############################################################################
#' @aliases show,taxonomyTable-method
#' @rdname show-methods
setMethod("show", "taxonomyTable", function(object){
    cat(paste("Taxonomy Table:     [", dim(object)[1], " taxa by ", 
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
	# print otu_table (always there).
	cat(paste("otu_table()   OTU Table:         [ ", ntaxa(otu_table(object)), " taxa and ", 
        nsamples(otu_table(object)), " samples ]", sep = ""), fill = TRUE)	

	# print Sample Data if there
	if(!is.null(sample_data(object, FALSE))){
        cat(paste("sample_data() Sample Data:       [ ", dim(sample_data(object))[1], " samples by ", 
	        dim(sample_data(object))[2], 
            " sample variables ]", sep = ""), fill = TRUE)
	}

	# print tax Tab if there	
	if(!is.null(tax_table(object, FALSE))){
        cat(paste("tax_table()   Taxonomy Table:    [ ", dim(tax_table(object))[1], " taxa by ", 
	        dim(tax_table(object))[2], 
            " taxonomic ranks ]", sep = ""), fill = TRUE)
	}
	
	# print tree if there
	if(!is.null(phy_tree(object, FALSE))){
        cat(paste("phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object)), " tips and ", 
	        phy_tree(object)$Nnode,
            " internal nodes ]", sep = ""),
        	fill = TRUE
        ) 
	}
	
	# print refseq summary if there
	if(!is.null(refseq(object, FALSE))){
        cat(paste("refseq()      ", class(refseq(object))[1], ":      [ ", ntaxa(refseq(object)), " reference sequences ]", sep = ""), fill=TRUE)
	}
	
})
############################################################################
