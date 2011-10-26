##############################################################################
# Extension of the transpose function in R for phyloseq objects.
# This is really only applied to the otuTable category.	
#
setGeneric("t")
setMethod("t", signature("otuTable"), function(x){
	#new("otuTable", t(x@.Data), speciesAreRows = (!speciesAreRows(x)))
	x <- otuTable( t(as(x, "matrix")), speciesAreRows=(!speciesAreRows(x)) )
	return(x)
})
##############################################################################
setMethod("t", signature("phyloseqFather"), function(x){
	x@otuTable <- t( otuTable(x) )
	return(x)
})
##############################################################################


