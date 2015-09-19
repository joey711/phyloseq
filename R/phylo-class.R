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
# Define horizontal position / node-ages by depth to root
# For instance, `xx` in `plot_tree` and `tipAges` in `fastUniFrac`
#' @keywords internal
ape_node_depth_edge_length <- function(Ntip, Nnode, edge, Nedge, edge.length){
  .C(ape:::node_depth_edgelength, PACKAGE="ape", as.integer(Ntip),
     as.integer(Nnode), as.integer(edge[, 1]),
     as.integer(edge[, 2]), as.integer(Nedge),
     as.double(edge.length), double(Ntip + Nnode))[[7]]
}