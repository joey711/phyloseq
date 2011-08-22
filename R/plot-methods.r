#
# extension of plot methods for different classes of objects in phyloseq package.
# 
################################################################################
#' Plot an otuTree object as a simple phylogenetic tree.
#'
#' @name plot
#' @aliases plot,otuTree-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "otuTree", function(x, ...){
	tree = tre(x)
	plot(tree, ...)	
	nodelabels(as.character(1:max(tree$edge)),node=1:max(tree$edge))
	edgelabels(as.character(1:nrow(tree$edge)),edge=1:nrow(tree$edge))
})
################################################################################
#' Plot a phyloseq object as a taxaplot style barplot.
#'
#' @name plot
#' @aliases plot,phyloseq-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "phyloseq", function(x, ...){
	taxaplot(otu=x, ...)
})
################################################################################
#' Plot a phyloseqTree object as an annotated phylogenetic
#'
#' @name plot
#' @aliases plot,phyloseqTree-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "phyloseqTree", function(x, ...){
	plot_tree_phyloseq(x, ...)
})
################################################################################
