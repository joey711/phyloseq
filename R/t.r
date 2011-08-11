##############################################################################
# Extension of the transpose function in R.
#
setGeneric("t")
setMethod("t", signature("otuTable"), function(x){
	new("otuTable", t(x@.Data),
		speciesAreRows=!speciesarerows(x))
})
##############################################################################





