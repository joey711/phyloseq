#
# extension of plot methods for phyloseq object.
# 
################################################################################
################################################################################
################################################################################
################################################################################
#' @export
#' @name plot
#' @aliases plot,phyloseq,ANY-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "phyloseq", function(x, ...){
	if( all(c("otuTable", "sampleMap", "tre") %in% names(splat.phyloseq.objects(x))) ){
		plot_tree_phyloseq(x, ...)		
	} else if( all(c("otuTable", "sampleMap", "taxTab") %in% names(splat.phyloseq.objects(x))) ){
		taxaplot(otu=x, ...)
	} else if( all(c("otuTable", "tre") %in% names(splat.phyloseq.objects(x))) ){
		tree <- as(tre(x), "phylo")
		ape::plot.phylo(tree, ...)	
		nodelabels(as.character(1:max(tree$edge)),node=1:max(tree$edge))
		edgelabels(as.character(1:nrow(tree$edge)),edge=1:nrow(tree$edge))		
	}
})
################################################################################