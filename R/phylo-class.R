# Methods related to using phylo in phyloseq, including 
# phyloseq-internal calls to ape internals.
################################################################################
#' Method for fixing problems with phylo-class trees in phyloseq
#' 
#' For now this only entails replacing each missing (\code{NA}) branch-length
#' value with 0.0.
#' 
#' @keywords internal
setGeneric("fix_phylo", function(tree) standardGeneric("fix_phylo") )
#' @rdname fix_phylo
#' @aliases fix_phylo,phylo-method
setMethod("fix_phylo", "phylo", function(tree){
  tree$edge.length[which(is.na(tree$edge.length))] <- 0
  return(tree)
})
################################################################################
