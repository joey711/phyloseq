#
# extension of plot methods for different classes of objects in phyloseq package.
# 
################################################################################
# # #' Plot an otuTree object as a simple phylogenetic tree.
# # #'
# # #' @param x an \code{otuTree} object.
# # #' @param ... Additional plotting parameters, providede to \code{ape::plot.phylo()}
#' @export
#' @name plot
#' @aliases plot,otuTree,ANY-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "otuTree", function(x, ...){
	tree <- as(tre(x), "phylo")
	ape::plot.phylo(tree, ...)	
	nodelabels(as.character(1:max(tree$edge)),node=1:max(tree$edge))
	edgelabels(as.character(1:nrow(tree$edge)),edge=1:nrow(tree$edge))
})
################################################################################
#' Plot a otuSam object as a taxaplot style barplot.
#'
#' @export
#' @name plot
#' @aliases plot,otuSamTax,ANY-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "otuSamTax", function(x, ...){
	taxaplot(otu=x, ...)
})
################################################################################
#' Plot a otuSamTree object as an annotated phylogenetic
#'
#' @export
#' @name plot
#' @aliases plot,otuSamTree,ANY-method
#' @docType methods
#' @rdname plot-methods
setMethod("plot", "otuSamTree", function(x, ...){
	plot_tree_phyloseq(x, ...)
})
################################################################################
