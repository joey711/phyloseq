################################################################################
# coercion methods
################################################################################
setAs("otuTable", "matrix", function(from){
	from@.Data
})
setAs("otuTree", "otuTable", function(from){
	otuTable(from)	
})
setAs("otuSam", "otuTable", function(from){
	otuTable(from)	
})
################################################################################
setAs("sampleMap", "data.frame", function(from){
	new("sampleMap", from)
})
setAs("data.frame", "sampleMap", function(from){
	data.frame(from)
})
setAs("otuSam", "sampleMap", function(from){
	sampleMap(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("otuSamTax", "taxonomyTable", function(from){
	taxTab(from)	
})
################################################################################
