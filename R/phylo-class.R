# phyloseq-specific definition of "phylo" class,
# and methods related to using phylo in phyloseq, including 
# phyloseq-internal calls to ape internals.
################################################################################
#
#' S3 class placeholder definition (list) for phylogenetic trees.
#' 
#' The ape package does not export a version of its \code{\link[ape]{phylo}}-class,
#' partly because it is not really defined formally anywhere.
#' Instead, it is an S3 class extended from the base class, \code{\link{list}} --
#' this is a very common and easy approach --
#' and proper behavior of any method taking an instance of this class 
#' requires exact naming conventions for element names of the components.
#' The phyloseq package does not provide any validity checks that a given phylo
#' instance is valid (conforms to the conventions in the ape package). Yet.
#' If problems arise, this might be considered, and they could be defined
#' judiciously and within phyloseq. 
#' Similarly, if a formal definition for the the phylo-class is ever exported
#' by ape, the current philosophy of phyloseq would be to remove this
#' internal definition and import the former. Note that there is still some 
#' work going on for the phylobase package, which is addressing these same 
#' exact issues for S4 phylogenetic tree interaction. 
#' A very large number of packages (around 60 at my last count), depend on ape,
#' making it easily the de facto standard for representing phylogenetic trees in R;
#' and the phyloseq team would prefer to use any exported definitions from
#' the ape package if possible and available.
#' 
#' @seealso 
#' \code{\link[ape]{phylo}}
#' 
#' @keywords internal
phylo <- structure(list(), class = "phylo")
################################################################################
# If this ever works
# @importClassesFrom ape phylo
################################################################################
#' An S4 placeholder of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree. It is used throught ape.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
#' @name phylo-class
#' @rdname phylo-class
#' @exportClass phylo
setOldClass("phylo")
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