# show.r
# 
# methods for displaying summaries of data on screen.
#
# 
############################################################################	
# print/show functions
# General specialized table-print function
############################################################################	
# nice looking dimnames summary function
dimsum = function(z){
	if( length(z) < 6 ){ 
		return(z)
	} else {
		paste(paste( head(z, min(c(3,length(z))) ),sep="",collapse=", "),
			" ... ",
			  paste( tail(z, min(c(3,length(z))) ),sep="",collapse=", "),
		sep="",collapse="")
	}
}
############################################################################	
tableprint = function(y, 
	typelabel="table", dim1label="rownames", dim2label="colnames"){
	### Print OTU table features.
	if( any(dim(y))){
		cat( paste(typelabel," [",dim(y)[1]," by ",dim(y)[2],"]:",sep=""),fill=TRUE)
		cat( paste(dim1label,":", sep=""), fill=TRUE)
		cat( dimsum(rownames(y)), fill=TRUE ) 
		cat( paste(dim2label,":",sep=""), fill=TRUE) 
		cat( dimsum(colnames(y)), fill=TRUE ) 
		# print tiny matrix part
		if( nrow(y) < 6 ){
			print( y[,1:min(c(5,ncol(y))),drop=FALSE] )
		} else {
			print( head(y[,1:min(c(5,ncol(y))),drop=FALSE], min(c(3,nrow(y)))) )
			cat("...",fill=TRUE)
			print( tail(y[,1:min(c(5,ncol(y))),drop=FALSE], min(c(3,nrow(y)))) )
		}
	}			
}
setMethod("print", "otuTable", function(x){
	# Prevent infinite subsetting loops by specifying .Data
	y = x@.Data		
	if(x@speciesAreRows){
		tableprint(y,"OTU Table","Species","Samples")			
	}else{
		tableprint(y,"OTU Table","Samples","Species")
	}			 
})
setMethod("show", "otuTable", function(object){print(object)})
setMethod("print", "sampleMap", function(x){
	tableprint(data.frame(x),"Sample Map","Samples","Variables")
})
setMethod("show", "sampleMap", function(object){print(object)})
setMethod("print", "taxonomyTable", function(x){
	y = x@.Data			
	tableprint(y,"Taxonomy Table","Species","Taxonomic Level")		
})
setMethod("show", "taxonomyTable", function(object){print(object)})

# For object printing at top of each show method for HO objects:
objectprint = function(x){
	cat(paste(class(x),"Object",sep=" "),fill=TRUE)
}
# print/show methods for the heterogenous objects.
setMethod("print", "phyloseq", function(x){
	objectprint(x)
	print(x@otuTable)
	cat("   ",fill=TRUE)
	print(x@sampleMap)
})
setMethod("show", "phyloseq", function(object){print(object)})

setMethod("print", "otuTree", function(x){
	objectprint(x)
	print(x@otuTable)
	cat("   ",fill=TRUE)
	print(x@tre)
})
setMethod("show", "otuTree", function(object){print(object)})

setMethod("print", "otuTree4", function(x){
  objectprint(x)
	print(x@otuTable)
	cat("   ",fill=TRUE)
	print(summary(x@tre))
})
setMethod("show", "otuTree4", function(object){print(object)})

setMethod("print", "phyloseqTax", function(x){
	objectprint(x)
	print(phyloseq(x))
	cat("   ",fill=TRUE)
	print(x@taxTab)
})
setMethod("show", "phyloseqTax", function(object){print(object)})

setMethod("print", "phyloseqTaxTree", function(x){
	objectprint(x)
	print(phyloseqTax(x))
	cat("   ",fill=TRUE)
	print(x@tre)
})
setMethod("show", "phyloseqTaxTree", function(object){print(object)})

setMethod("print", "phyloseqTree", function(x){
	objectprint(x)
	print(phyloseq(x))
	cat("   ",fill=TRUE)
	print(x@tre)
})
setMethod("show", "phyloseqTree", function(object){print(object)})
############################################################################	