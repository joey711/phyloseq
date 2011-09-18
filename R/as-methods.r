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
setAs("data.frame", "sampleMap", function(from){
	new("sampleMap", from)
})
setAs("sampleMap", "data.frame", function(from){
	data.frame(from)
})
setAs("otuSam", "sampleMap", function(from){
	sampleMap(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("otuTax", "taxonomyTable", function(from){
	taxTab(from)	
})
################################################################################
setAs("otuTree", "phylo4", function(from){
	tre(from)	
})
setAs("otuTree", "phylo", function(from){
	as(tre(from), "phylo")
})
################################################################################