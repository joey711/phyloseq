################################################################################
#' Build objects for the phyloseq package.
#'
#' \code{phyloseq()} is both a constructor and accessor method, on account of
#' the ambiguation of `phyloseq' being both the package name and a class
#' name (the \code{phyloseq} class, which contains an \code{otuTable} and 
#' \code{sampleMap} slots, only).
#'
#' This is the main method suggested for constructing \code{phyloseq} higher-order
#' objects from their components.
#'
#' @param object An object among the set of classes defined by the phyloseq package,
#' as well as \code{phylo} and \code{phylo4} objects (which are defined by the
#' ape and phylobase packages, respectively). For constructing H.O. objects,
#' this should be an otuTable object.
#' @param x1 (optional) an additional component data object (but not otuTable)
#' @param x2 (optional) an additional component data object (but not otuTable)
#' @param x3 (optional) an additional component data object (but not otuTable)
#'
#' @return The class of the returned object depends on the argument 
#' class(es). To construct a H.O. object, the otuTable must be the first object
#' and/or assigned to the named argument \code{object}.
#' Otherwise, the order of arguments does not matter. If a single component-class
#' object is provided, it is simply returned as-is. If a higher-order object is
#' provided as input, then the accessor method is invoked, and an object of 
#' class \code{phyloseq} is returned.
#'
#' @seealso merge_phyloseq
#' @keywords phyloseq constructor
#' @export
#' @examples #
#' ## OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
#' ## tax1 <- taxTab(matrix("abc", 30, 8))
#' ## map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' ##   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) ) 
#' ## map1 <- sampleMap(map1)
#' ## ex1 <- phyloseq(OTU1, map1, tax1)
#' ## phyloseq(ex1)
setGeneric("phyloseq", function(object, x1, x2, x3) standardGeneric("phyloseq"))
setMethod("phyloseq", "phyloseq", function(object){
	new("phyloseq", otuTable=otuTable(object), sampleMap=sampleMap(object))		
})
setMethod("phyloseq", "phyloseqTax", function(object){ callNextMethod(object) })
setMethod("phyloseq", "phyloseqTree", function(object){ callNextMethod(object) })
setMethod("phyloseq", "phyloseqTaxTree", function(object){ callNextMethod(object) })
################################################################################
# The set of methods where the argument is a single component data object: 
# just return object.
setMethod("phyloseq", "otuTable",      function(object){ return(object) })
setMethod("phyloseq", "sampleMap",     function(object){ return(object) })
setMethod("phyloseq", "taxonomyTable", function(object){ return(object) })
setMethod("phyloseq", "phylo",         function(object){ return(object) })
setMethod("phyloseq", "phylo4",        function(object){ return(object) })
################################################################################
setMethod("phyloseq", signature("otuTable", "sampleMap"), function(object, x1){
  new("phyloseq", otuTable=object, sampleMap=x1)
})

# otuTree / otuTree4
setMethod("phyloseq", signature("otuTable", "phylo"), function(object, x1){
  new("otuTree", otuTable=object, tre=x1)
})
setMethod("phyloseq", signature("otuTable", "phylo4"), function(object, x1){
  new("otuTree4", otuTable=object, tre=x1)
})

# phyloseqTax
setMethod("phyloseq", signature("otuTable", "sampleMap", "taxonomyTable"),
          function(object, x1, x2){
  new("phyloseqTax", otuTable=object, sampleMap=x1, taxTab=x2)
})
setMethod("phyloseq", signature("otuTable", "taxonomyTable", "sampleMap"),
          function(object, x1, x2){
  new("phyloseqTax", otuTable=object, sampleMap=x2, taxTab=x1)
})
  
# phyloseqTree
setMethod("phyloseq", signature("otuTable", "phylo", "sampleMap"),
          function(object, x1, x2){
  new("phyloseqTree", otuTable=object, tre=x1, sampleMap=x2)
})
setMethod("phyloseq", signature("otuTable", "sampleMap", "phylo"),
          function(object, x1, x2){
  new("phyloseqTree", otuTable=object, tre=x2, sampleMap=x1)
})
  
# phyloseqTaxTree
setMethod("phyloseq", signature("otuTable", "sampleMap", "phylo", "taxonomyTable"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x1, tre=x2, taxTab=x3)
})
setMethod("phyloseq", signature("otuTable", "sampleMap", "taxonomyTable", "phylo"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x1, tre=x3, taxTab=x2)
})
setMethod("phyloseq", signature("otuTable", "phylo", "taxonomyTable", "sampleMap"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x3, tre=x1, taxTab=x2)
})
setMethod("phyloseq", signature("otuTable", "phylo", "sampleMap", "taxonomyTable"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x2, tre=x1, taxTab=x3)
})
setMethod("phyloseq", signature("otuTable", "taxonomyTable", "sampleMap", "phylo"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x2, tre=x3, taxTab=x1)
})
setMethod("phyloseq", signature("otuTable", "taxonomyTable", "phylo", "sampleMap"),
          function(object, x1, x2, x3){
  new("phyloseqTaxTree", otuTable=object, sampleMap=x3, tre=x2, taxTab=x1)
})
################################################################################
# tests:
#phyloseq(otuTable(ex4), sampleMap(ex4), taxTab(ex4), tre(ex4) )
#phyloseq(otuTable(ex4), sampleMap(ex4), taxTab(ex4) )
#phyloseq(otuTable(ex4), sampleMap(ex4), tre(ex4) )
#phyloseq(otuTable(ex4), tre(ex4) )
#phyloseq(otuTable(ex4), sampleMap(ex4) )