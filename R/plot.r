#
# extension of plot methods for different classes of objects in phyloseq package.
# 
################################################################################
setMethod("plot", "otuTree", function(x, ...){
	#require(ape)
	tree = tre(x)
	plot(tree, ...)	
	nodelabels(as.character(1:max(tree$edge)),node=1:max(tree$edge))
	edgelabels(as.character(1:nrow(tree$edge)),edge=1:nrow(tree$edge))
})
################################################################################
setMethod("plot", "phyloseq", function(x, ...){
	# Should be a nice default of phylabarplot
  taxaplot(otu=x, ...)
  #list(x, ...)
})
################################################################################


