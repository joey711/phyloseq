################################################################################
#' Agglomerate closely-related taxa using single-linkage clustering.
#' 
#' All tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated into one taxa using \code{merge_species}.
#' 
#' Can be used to create a non-trivial OTU Table, if a phylogenetic tree is available.
#'
#' For now, a simple, ``greedy'', single-linkage clustering is used. In future releases
#' it should be possible to specify different clustering approaches available in \code{R},
#' in particular, complete-linkage clustering appears to be used more commonly for OTU
#' clustering applications.
#'
#' @usage tipglom(tree, OTU, speciationMinLength=0.02)
#'
#' @param tree \code{\link{phyloseq-class}}, containing an OTU Table and
#'  phylogenetic tree. If, alternatively, \code{tree} is a \code{\link{phylo-class}},
#'  then \code{OTU} is required.
#'
#' @param OTU An otuTable object. Optional. Ignored if \code{tree} is a 
#'  \code{\link{phyloseq-class}} object. If \code{tree} is a \code{phylo}
#'  object and \code{OTU} is provided, then return will be an \code{phyloseq}
#'  object. 
#'
#' @param speciationMinLength The minimum branch length that separates taxa. All
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated. Default is 0.02
#'
#' @return An object of class \code{phyloseq}. Output class matches
#' the class of \code{tree}, unless it is a \code{phylo} object, in
#' which case \code{tipglom} returns an \code{phyloseq} object.
#'
#' @rdname tipglom-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # # data(phylocom)
#' # # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # # # x1   <- phyloseq(otu, phylocom$phylo)
#' # # # print(x1); par(mfrow=c(2, 1)); plot(tre(x1))
#' # # # x2 <- tipglom(x1, speciationMinLength = 2.5)
#' # # # plot(tre(x2))
#' # # # ## Try on example datset 1
#' # # # data(ex1); nspecies(ex1)
#' # # # ex7 <- tipglom(ex1, speciationMinLength = 0.05)
#' # # # nspecies(ex7)
#' # data(esophagus); nspecies(esophagus); par(mfrow=c(2, 1)); plot(tre(esophagus))
#' # tre(esophagus)$edge.length
#' # x3 <- tipglom(esophagus, speciationMinLength = 0.20)
#' # nspecies(x3); plot(tre(x3))
setGeneric("tipglom", function(tree, OTU, speciationMinLength=0.02) standardGeneric("tipglom"))
#' @rdname tipglom-methods
#' @aliases tipglom,phylo,otuTable-method
setMethod("tipglom", signature("phylo", "otuTable"), function(tree, OTU, speciationMinLength=0.02){
	# dispatch as the combined-object (phyloseq-class), auto coherence.
	tipglom( phyloseq(OTU, tree), speciationMinLength=speciationMinLength)
})
#' @rdname tipglom-methods
#' @aliases tipglom,phyloseq,ANY-method
setMethod("tipglom", signature("phyloseq"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
#' @rdname tipglom-methods
#' @aliases tipglom,phylo,ANY-method
setMethod("tipglom", signature("phylo"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
################################################################################
#' Internal function for tiplgom.
#' 
#' Internal function, users should use the S4 method \code{\link{tipglom}}.
#' Tree can be a \code{\link{phyloseq-class}} that contains a phylogenetic tree, 
#' This is because \code{\link{merge_species}} can
#' handle all the relevant objects, as can \code{\link{getTipDistMatrix}}.
#' Create Non-trivial OTU table, by agglomerating nearby tips.
#' tipglom.internal is called by the S4 \code{tipglom} methods. It is useful if 
#' a motivated user wants to see the internals of the implementation. By design
#' it lacks explicit object handling. Use \code{\link{tipglom}} instead.
#'
#' @param tree An object of class \code{phylo}, or \code{phyloseq} 
#'
#' @param speciationMinLength The minimum branch length that separates taxa. All
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated.
#'
#' @return An object of class \code{phylo}, or \code{phyloseq}.
#'  Output class matches the class of \code{tree}.
#'
#' @seealso tipglom
#' @import igraph
#' @keywords internal
tipglom.internal <- function(tree, speciationMinLength){
	# Create adjacency matrix, where tips are adjacent
	# if their distance is below the threshold speciationMinLength
	tipAdjacent <- (getTipDistMatrix( tree ) < speciationMinLength)
	# Define igraph object based on the tree-tip adjacenecy matrix
	ig          <- graph.adjacency(tipAdjacent, diag=FALSE)
	# Define the species cliques to loop through
	spCliques   <- edgelist2clique( get.edgelist(ig) )
	# successively merge
	for( i in 1:length(spCliques)){
		tree <- merge_species(tree, eqspecies=spCliques[[i]])
	}
	# Test if you missed anything:
	# graph.adjacency( getTipDistMatrix(tree) < speciationMinLength, diag=FALSE )
	return(tree)
}
#################################################################
#' An internal wrapper function on \code{\link[ape]{cophenetic.phylo}}
#' 
#' This is useful for determining tips that should be combined.
#' 
#' @param tree \code{phylo}
#' 
#' @param byRootFraction Should the distance be calculated according to
#' fractional distance to the root? If \code{FALSE}, the distance is
#' instead the patristic distance as calculated by cophenetic.phylo. 
#' Default \code{FALSE}.
#' 
#' @return character matrix. First column is the complete match, followed by
#'   one for each capture group
#' 
#' @seealso tipglom
#' @keywords internal
#' @aliases gettipdistmatrix getTipDistMatrix
setGeneric("getTipDistMatrix", function(tree, byRootFraction=FALSE) standardGeneric("getTipDistMatrix"))
setMethod("getTipDistMatrix", signature("phylo"), function(tree, byRootFraction=FALSE){
	### require("picante") # picante is a "depends"-level dependency of phyloseq.
	pairwiseSpecDists = cophenetic(tree)
	# If byRootFraction is true, normalize the cophenetic distances
	# according to the mean root age.
	if( byRootFraction ){
		# Want to normalize pairwise tip distances by their mean distance to root
		# start with tipAges
		tipAges = node.age(tree)$ages[which(tree$edge[,2] %in% 1:length(tree$tip.label))]
		names(tipAges) = tree$tip.label
		###### Want Mmean to be a matrix of the mean pairwise root-distance b/w each tip-pair
		Mmean = matrix(NA,length(tipAges),length(tipAges),
			dimnames=list(names(tipAges),names(tipAges))) 	
		means = combn(tipAges,2,mean)
		ind = combn(length(tipAges),2)
		for(i in 1:ncol(ind)){Mmean[ind[1,i], ind[2,i]] <- means[i]}
		for(i in 1:ncol(ind)){Mmean[ind[2,i], ind[1,i]] <- means[i]}
		diag(Mmean) <- tipAges
		# take the ratio of spec distances to the mean
		fracDists = pairwiseSpecDists / Mmean
		return(fracDists)
	} else {
		return(pairwiseSpecDists)
	}
})
setMethod("getTipDistMatrix", signature("phyloseq"), function(tree, byRootFraction=FALSE){
	getTipDistMatrix( tre(tree) )
})
gettipdistmatrix <- getTipDistMatrix
################################################################################
#' Convert edgelist hash-table to clique list
#'
#' Agglomerating function to convert an edgelist -- which is a 2-column table
#' of vertices where each row represents an edge -- to a list of cliques,
#' in which each clique is represented by a character vector of the vertex labels
#' of the vertices that are members of the clique. This algorithm is perfectly
#' greedy, such that the only requirement for inclusion in a clique is an edge
#' to any of the other members of that clique.
#'
#' @usage edgelist2clique(EdgeList)
#'
#' @param EdgeList a 2-column table of vertices where each row represents an edge. 
#'
#' @return A list, where each element is a character vector of tips that should
#' are in the same clique.
#'
#' @keywords internal
#' @examples #
#' # edgelist2clique(get.edgelist(ig))
edgelist2clique = function(EdgeList){
	# initialize list of globs
	glob	= vector(mode="list")
	
	for (i in 1:nrow(EdgeList) ){
		# initialize for each loop a 'skip' variable, to avoid later tests if answer already found.
		skip	= FALSE
		# check entries in list, glob, for membership in already-forming glob
		thislink	= unlist(EdgeList[i,1:2])
		# identify which, if any, globs have contig of interest already in them
		# OLD WAY 1
		# glob1	= which(sapply(sapply(glob,function(globi,ctig){which(ctig==globi)},ctig=thislink[1]),length)>0)
		# glob2	= which(sapply(sapply(glob,function(globi,ctig){which(ctig==globi)},ctig=thislink[2]),length)>0)
		# better way
		glob1 = which(vapply(glob, function(glb,vertex){vertex %in% glb}, TRUE, thislink[1]))
		glob2 = which(vapply(glob, function(glb,vertex){vertex %in% glb}, TRUE, thislink[2]))
		# grep is actually slower
		# glob1 = grep(thislink[1], glob, fixed=TRUE)
		# glob2 = grep(thislink[2], glob, fixed=TRUE)
	
		## Now series of if-tests to decide where and how to glob.
		# if both contigs are already in same globs, skip
		if( !skip & length(glob1)>0 & length(glob2)>0 ){
			if( glob1==glob2 ){
				# skip remaining tests, only.
				skip = TRUE	
			}
		}
		# if both contigs are already in different globs, join both globs
		if( !skip & length(glob1)>0 & length(glob2)>0 ){
			if( glob1!=glob2 ){
				# join globs at end	
				glob	= c(glob,list(c(glob[[glob1]],glob[[glob2]])))
				# remove old globs
				glob	= glob[-c(glob1,glob2)]
				# skip remaining tests.
				skip = TRUE	
			}
		}
		# if only contig 1 is in a glob, add Contig 2 to that glob
		if( !skip & length(glob1)>0 & length(glob2)==0){
			# add Contig 2 to glob1	
			glob[[glob1]]	= c(glob[[glob1]],thislink[2])
			# skip remaining tests.
			skip = TRUE	
		}
		# if only Contig 2 is in a glob, add contig 1 to that glob
		if( !skip & length(glob2)>0 & length(glob1)==0 ){
			# add Contig 1 to glob2	
			glob[[glob2]]	= c(glob[[glob2]],thislink[1])
			# skip remaining tests.
			skip = TRUE		
		}
		# all else, form new glob
		if( !skip & length(glob1)==0 & length(glob2)==0 ){
			# add Contig 1 and Contig 2 as vector in new glob.
			glob	= c(glob,list(thislink))
		}
	}
	return(glob)
}
################################################################################
################################################################################
################################################################################
################################################################################
#' Agglomerate taxa of the same type.
#'
#' This method merges species if, at a certain taxaonomic rank, their taxonomy
#' is the same. 
#' Its approach is analogous to \code{tipglom}, but uses categorical data
#' instead of a tree. In principal, other categorical data known for all taxa
#' could also be used in place of taxonomy. 
#'
#' @usage taxglom(physeq, tax=NULL, taxlevel="Phylum", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} or \code{\link{otuTable}}.
#'
#' @param tax (Optional). Either a \code{link{taxonomyTable-class}}, or alternatively, a 
#'  character vector specifying the desired taxonomic group of each taxa in 
#'  \code{physeq}. If \code{tax} is a character vector, it must have length equal
#'  to the (original) number of taxa in \code{physeq} (\code{nspecies(physeq)}), 
#'  and each element must be
#'  named according to the taxa ID (that is, the result of 
#'  \code{species.names(physeq)}). If \code{tax} is a character vector, than
#'  the \code{taxlevel} argument is ignored. If \code{physeq} already contains
#'  a \code{taxonomyTable} component in its \code{taxTab} slot, then 
#'  the \code{tax} argument is ignored. 
#'
#' @param taxlevel A single-element character specifying the taxonomic level
#'  (column name)
#'  in \code{tax}, the \code{taxonomyTable}, that you want to agglomerate over.
#'  The default value is \code{"Phylum"}. Note that this default may
#'  agglomerate too broadly for a given experiment, and the user is strongly
#'  encouraged to try different taxonomic levels.
#'
#' @param NArm (Optional). Logical, length equal to one. Default is \code{TRUE}.
#'  CAUTION. The decision to prune (or not) taxa for which you lack categorical
#'  data could have a large effect on downstream analysis. You may want to
#'  re-compute your analysis under both conditions, or at least think carefully
#'  about what the effect might be and the reasons explaining the absence of 
#'  information for certain taxa. In the case of taxonomy, it is often a result 
#'  of imprecision in taxonomic designation based on short phylogenetic sequences
#'  and a patchy system of nomenclature. If this seems to be an issue for your
#'  analysis, think about also trying the nomenclature-agnostic \code{\link{tipglom}}
#'  method if you have a phylogenetic tree available.
#'
#' @param bad_empty (Optional). Character vector. Default: \code{c(NA, "", " ", "\t")}.
#'  Defines the bad/empty values 
#'  that should be ignored and/or considered unknown. They will be removed
#'  from the internal agglomeration vector derived from the argument to \code{tax},
#'  and therefore agglomeration will not combine taxa according to the presence
#'  of these values in \code{tax}. Furthermore, the corresponding taxa can be
#'  optionally pruned from the output if \code{NArm} is set to \code{TRUE}.
#' 
#' @return A taxonomically-agglomerated, optionally-pruned, object with class matching
#' the class of \code{physeq}.
#'
#' @seealso \code{\link{tipglom}}, \code{\link{prune_species}}, \code{\link{merge_species}}
#' 
#' @rdname taxglom-methods
#' @docType methods
#' @export
#'
#' @examples
#' # data(ex1)
#' # ## print the available taxonomic ranks
#' # colnames(taxTab(ex1))
#' # ## agglomerate at the Family taxonomic rank
#' # (x1 <- taxglom(ex1, taxlevel="Family") )
#' # ## How many taxa before/after agglomeration?
#' # nspecies(ex1); nspecies(x1)
#' # ## Look at enterotype dataset...
#' # data(enterotype)
#' # ## print the available taxonomic ranks. Shows only 1 rank available, not useful for taxglom
#' # colnames(taxTab(enterotype))
setGeneric("taxglom", 
	function(physeq, tax=NULL, taxlevel="Phylum", NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){
	standardGeneric("taxglom")
})
#' @rdname taxglom-methods
#' @aliases taxglom,otuTable,taxonomyTable-method
setMethod("taxglom", c("otuTable", "taxonomyTable"),
				function(physeq, tax=NULL, taxlevel="Phylum", NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){
	# vectorize the taxonomy table.
	tax <- as(tax, "matrix")[, taxlevel]
	taxglom.internal(physeq, tax, NArm, bad_empty)
})
#' @rdname taxglom-methods
#' @aliases taxglom,otuTable,character-method
setMethod("taxglom", c("otuTable", "character"),
				function(physeq, tax=NULL,taxlevel="Phylum", NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){
	taxglom.internal(physeq, tax, NArm, bad_empty)
})
#' @rdname taxglom-methods
#' @aliases taxglom,phyloseq,ANY-method
setMethod("taxglom", "phyloseq",
				function(physeq, tax=NULL, taxlevel="Phylum", NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){
	# vectorize the taxonomy table.
	tax <- as(taxTab(physeq), "matrix")[, taxlevel]
	taxglom.internal(physeq, tax, NArm, bad_empty)
})
################################################################################
#' taxglom core internal function.
#'
#' taxglom.internal makes all the glomming happen, and delegates the
#' object-handling issues to \code{merge_species()}.
#'
#' @param physeq the object on which agglomeration is to take place.
#'
#' @param tax for this internal function, tax must be a character vector.
#' \code{tax} is a vector in the core internal function.
#' 
#' @return the agglomerated object. Class matches argument \code{physeq}.
#'
#' @keywords internal
taxglom.internal <- function(physeq, tax, NArm=TRUE, bad_empty=c(NA, "", " ", "\t") ){

	# if NArm is TRUE, remove the empty, white-space, NA values from 
	if( NArm ){
		keep_species <- names(tax)[ !(tax %in% bad_empty) ]
		physeq <- prune_species(keep_species, physeq)
	}

	# Remove NAs and useless from the vector/factor for looping.
	# This does not remove the taxa that have an unknown (NA)
	# taxonomic designation at this particular taxonomic rank.
	tax <- tax[ !(tax %in% bad_empty) ]
	
	# Define the species cliques to loop through
	spCliques <- tapply(names(tax), factor(tax), list)
	
	# successively merge taxa in physeq.
	for( i in names(spCliques)){
		# print(i)
		physeq <- merge_species(physeq, eqspecies=spCliques[[i]])
	}
	return(physeq)
}
################################################################################
# test <- taxglom.internal(ex1, as(taxTab(ex1), "matrix")[, "Phylum"])
# testvec <- as(taxTab(ex1), "matrix")[, "Phylum", drop=TRUE]
# tapply(names(testvec), factor(testvec), length)
################################################################################
#' Prune unwanted species / taxa from a phylogenetic object.
#' 
#' An S4 Generic method for removing (pruning) unwanted taxa from phylogenetic
#' objects, including phylo-class trees, as well as native phyloseq package
#' objects. This is particularly useful for pruning a phyloseq object that has
#' more than one component that describes species.
#' The \code{phylo}-class version is adapted from \code{picante::prune.samples}.
#'
#' @param species (Required). A character vector of the species in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept species are TRUE, and length
#' is equal to the number of species in object x. If \code{species} is a named
#' logical, the species retained is based on those names. Make sure they are
#' compatible with the \code{species.names} of the object you are modifying (\code{x}). 
#'
#' @param x (Required). A phylogenetic object, including \code{phylo} trees,
#' as well as all phyloseq classes that represent taxa / species. If the function
#' \code{\link{species.names}} returns a non-\code{NULL} value, then your object
#' can be pruned by this function.
#'
#' @return The class of the object returned by \code{prune_species} matches
#' the class of the argument, \code{x}.
#'
#' @name prune_species
#' @rdname prune_species-methods
#' @export
#' @examples #
#' ## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' ## f1  <- filterfunSample(topk(2))
#' ## wh1 <- genefilterSample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_species(wh1, testOTU)
#' ## prune_species(wh2, testOTU)
#' ## 
#' ## taxtab1 <- taxTab(matrix("abc", 5, 5))
#' ## prune_species(wh1, taxtab1)
#' ## prune_species(wh2, taxtab1)
setGeneric("prune_species", function(species, x) standardGeneric("prune_species"))
################################################################################
#' @aliases prune_species,NULL,ANY-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("NULL"), function(species, x){
	return(x)
})
################################################################################
#' @aliases prune_species,character,phylo-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("character", "phylo"), function(species, x){
	trimTaxa <- setdiff(x$tip.label, species)
	if( length(trimTaxa) > 0 ){
		ape::drop.tip(x, trimTaxa)
	} else x
})
################################################################################
#' @aliases prune_species,character,otuTable-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("character", "otuTable"), function(species, x){
	species <- intersect( species, species.names(x) )
	if( speciesarerows(x) ){
		x[species, , drop=FALSE]
	} else {
		x[, species, drop=FALSE]
	}	
})
################################################################################
#' @aliases prune_species,character,sampleData-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("character", "sampleData"), function(species, x){
	return(x)
})
################################################################################
#' @aliases prune_species,character,phyloseq-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("character", "phyloseq"), 
		function(species, x){
			
	# Save time and return if the union of all component species names
	# captured by species.names(x) is same as species. 
	if( setequal(species.names(x), species) ){
		return(x)
	} else {	
		# All phyloseq objects have an otuTable slot, no need to test.
		x@otuTable   <- prune_species(species, otuTable(x))
		
		# Test if slot is present. If so, perform the component prune.
		if( !is.null(access(x, "taxTab")) ){
			x@taxTab <- prune_species(species, taxTab(x))
		}
		if( !is.null(access(x, "tre")) ){
			x@tre    <- prune_species(species, tre(x))
		}
		return(x)
	}
})
################################################################################
#' @aliases prune_species,character,taxonomyTable-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("character", "taxonomyTable"), 
		function(species, x){
	species <- intersect( species, species.names(x) )
	return( x[species, , drop=FALSE] )
})
################################################################################
#' @aliases prune_species,logical,ANY-method
#' @rdname prune_species-methods
setMethod("prune_species", signature("logical", "ANY"), function(species, x){
	# convert the logical argument to character and dispatch
	if( is.null(names(species)) ){
		species <- species.names(x)[species]
	} else {
		species <- names(species)[species]
	}
	prune_species(species, x)
})
################################################################################
################################################################################
#' Prune unwanted samples from a phyloseq object.
#' 
#' An S4 Generic method for removing (pruning) unwanted samples.
#'
#' @usage prune_samples(samples, x)
#'
#' @param samples A character vector of the samples in object x that you want to
#' keep. 
#'
#' @param x A phyloseq object.
#'
#' @return The class of the object returned by \code{prune_samples} matches
#' the class of the phyloseq object, \code{x}.
#'
#' @seealso \code{\link{subset_samples}}
#' 
#' @rdname prune_samples-methods
#' @docType methods
#' @export
#' @examples #
#'  # data(ex1)
#'  # B_only_sample_names <- sample.names(sampleData(ex1)[(sampleData(ex1)[, "Gender"]=="B"),])
#'  # ex2 <- prune_samples(B_only_sample_names, ex1)
#'  # ex3 <- subset_samples(ex1, Gender=="B")
#'  # ## This should be TRUE.
#'  # identical(ex2, ex3)
#'  # ## Here is a simpler example: Make new object with only the first 5 samples
#'  # ex4 <- prune_samples(sample.names(ex1)[1:5], ex1)
setGeneric("prune_samples", function(samples, x) standardGeneric("prune_samples"))
################################################################################
#' @aliases prune_samples,character,otuTable-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "otuTable"), function(samples, x){
	if( speciesarerows(x) ){
		x[, samples]
	} else {
		x[samples, ]
	}
})
################################################################################
#' @aliases prune_samples,character,sampleData-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "sampleData"), function(samples, x){
	x[samples, ]
})
################################################################################
#' @aliases prune_samples,character,phyloseq-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "phyloseq"), function(samples, x){
	x@samData  <- prune_samples(samples, sampleData(x) )
	x@otuTable <- prune_samples(samples, otuTable(x) )
	return(x)
})
################################################################################
####################################################################################
#' Thresholded rank transformation.
#' 
#' The lowest \code{thresh} values in \code{x} all get the value 'thresh'.
#'
#' @usage threshrank(x, thresh, keep0s=FALSE, ...)
#'
#' @param x (Required). Numeric vector to transform.
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain 
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#' 
#' @return A ranked, (optionally) thresholded numeric vector with length equal to
#'  \code{x}. Default arguments to \code{rank} are used, unless provided as
#'  additional arguments. 
#'
#' @seealso \code{\link{transformsamplecounts}}, \code{\link{rank}}, \code{\link{threshrankfun}}
#' @export 
#' @examples #
#' (a_vector <- sample(0:10, 100, TRUE))
#' threshrank(a_vector, 5, keep0s=TRUE)
#' data(ex1)
#' ## These three approaches result in identical otuTable
#' (x1 <- transformsamplecounts( otuTable(ex1), threshrankfun(500)) )
#' (x2 <- otuTable(apply(otuTable(ex1), 2, threshrankfun(500)), speciesAreRows(ex1)) )
#' identical(x1, x2)
#' (x3 <- otuTable(apply(otuTable(ex1), 2, threshrank, thresh=500), speciesAreRows(ex1)) )
#' identical(x1, x3)
threshrank <- function(x, thresh, keep0s=FALSE, ...){
	if( keep0s ){ index0 <- which(x == 0) }
	x <- rank(x, ...)
	thresh <- thresh[1]
	x[x<thresh] <- thresh
	if( keep0s ){ x[index0] <- 0 }
	return(x)
}
####################################################################################
#' A closure version of the \code{threshrank} function.
#'
#' Takes the same arguments as \code{\link{threshrank}}, except for \code{x}, 
#' because the output is a single-argument function rather than a rank-transformed numeric. 
#' This is useful for higher-order functions that require a single-argument function as input,
#' like \code{\link{transformsamplecounts}}.
#'
#' @usage threshrankfun(thresh, keep0s=FALSE, ...)
#' 
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain 
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#' 
#' @return A single-argument function with the options to \code{\link{threshrank}} set.
#'  
#' @seealso \code{\link{transformsamplecounts}}, \code{\link{threshrankfun}},
#'  \code{\link{threshrank}}
#' @export
#' @examples
#' data(ex1)
#' ## These three approaches result in identical otuTable
#' (x1 <- transformsamplecounts( otuTable(ex1), threshrankfun(500)) )
#' (x2 <- otuTable(apply(otuTable(ex1), 2, threshrankfun(500)), speciesAreRows(ex1)) )
#' identical(x1, x2)
#' (x3 <- otuTable(apply(otuTable(ex1), 2, threshrank, thresh=500), speciesAreRows(ex1)) )
#' identical(x1, x3)
threshrankfun <- function(thresh, keep0s=FALSE, ...){
	function(x){
		threshrank(x, thresh, keep0s=FALSE, ...)
	}
}
################################################################################
#' Transpose \code{\link{otuTable-class}} or \code{\link{phyloseq-class}}
#'
#' Extends the base transpose method, \code{\link[base]{t}}.
#'
#' @usage t(x)
#'
#' @param x An \code{otuTable} or \code{\link{phyloseq-class}}.
#'
#' @return The class of the object returned by \code{t} matches
#' the class of the argument, \code{x}. The \code{otuTable} is
#' transposed, and \code{\link{speciesAreRows}} value is toggled.
#'
#' @name t
#' @rdname transpose-methods
#' @docType methods
#' @export
#' @examples
#' data(ex1)
#' otuTable(ex1)
#' t( otuTable(ex1) )
setGeneric("t")
#' @aliases t,otuTable-method
#' @rdname transpose-methods
setMethod("t", signature("otuTable"), function(x){
	#new("otuTable", t(x@.Data), speciesAreRows = (!speciesAreRows(x)))
	x <- otuTable( t(as(x, "matrix")), speciesAreRows=(!speciesAreRows(x)) )
	return(x)
})
################################################################################
#' @aliases t,phyloseq-method
#' @rdname transpose-methods
setMethod("t", signature("phyloseq"), function(x){
	x@otuTable <- t( otuTable(x) )
	return(x)
})
################################################################################
#' Transform the abundance count data in an \code{otuTable}, sample-by-sample.
#' 
#' This function transforms the sample counts of a species
#' abundance matrix according to a user-provided function.
#' The counts of each sample will be transformed individually. No sample-sample 
#' interaction/comparison is possible by this method. 
#'
#' @usage transformsamplecounts(physeq, fun)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} of \code{\link{otuTable-class}}.
#'
#' @param fun (Required). A single-argument function that will be applied
#'  to the abundance counts of each sample. Can be an anonymous \code{\link[base]{function}}.
#' 
#' @return A transformed \code{otuTable} -- or \code{phyloseq} object with its
#'  transformed \code{otuTable}. 
#'  In general, trimming is not expected by this 
#'  method, so it is suggested that the user provide only functions that return
#'  a full-length vector. Filtering/trimming can follow, for which the 
#'  \code{\link{genefilterSample}} and \code{\link{prune_species}} functions
#'  are suggested.
#'
#' @seealso \code{\link{threshrankfun}}, \code{\link{rank}}, \code{\link{log}}
#'
#' @docType methods
#' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @rdname transformcounts
#' @export
#'
#' @examples #
#' data(ex1)
#' ## transformsamplecounts can work on phyloseq-class, modifying otuTable only
#' (ex1r <- transformsamplecounts(ex1, rank) )
#' ## These two approaches result in identical otuTable
#' (x1 <- transformsamplecounts( otuTable(ex1), threshrankfun(500)) )
#' (x2 <- otuTable(apply(otuTable(ex1), 2, threshrankfun(500)), speciesAreRows(ex1)) )
#' identical(x1, x2)
transformsamplecounts <- function(physeq, fun){
	if( speciesarerows(physeq) ){
		newphyseq <- apply(as(otuTable(physeq), "matrix"), 2, fun)
	} else {
		newphyseq <- apply(as(otuTable(physeq), "matrix"), 1, fun)
	}
	otuTable(physeq) <- otuTable(newphyseq, speciesAreRows=speciesarerows(physeq))
	return(physeq)
}
####################################################################################
# # #' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @docType methods
#' @rdname transformcounts
#' @export
TransformSampleCounts <- transformsamplecounts
####################################################################################
# # #' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @docType methods
#' @rdname transformcounts
#' @export
transformSampleCounts <- transformsamplecounts
####################################################################################
############################################################
#' Filter OTUs with arbitrary function, sample-wise.
#' 
#' A general OTU trimming function for selecting OTUs that satisfy
#' some criteria within the distribution of each sample, and then
#' also an additional criteria for number of samples that must pass.
#' This is a genefilter-like function that only considers sample-wise
#' criteria. The number of acceptable samples is used
#' as the final criteria (set by the argument \code{A})
#' to determine whether or not the taxa should
#' be retained (\code{TRUE}) or not (\code{FALSE}). Just like with genefilter, a 
#' logical having length equal to nrow()/\code{\link{nspecies}} is returned, indicating which
#' should be kept. This output can be provided
#' directly to OTU trimming function, \code{\link{prune_species}}.
#' By contrast, \code{\link[genefilter]{genefilter}}, 
#' of the genefilter package in Bioconductor,
#' works only on the rows of a matrix. Note that, because \code{\link{otuTable-class}}
#' inherits directly from the \code{\link{matrix-class}}, an unmodified
#' otuTable can be provided to \code{genefilter}, but be mindful of the orientation
#' of the otuTable (use \code{\link{speciesAreRows}}),
#' and transpose (\code{\link[phyloseq]{t}}) if needed.
#'
#' @usage genefilterSample(X, flist, A=1)
#'
#' @param X The object that needs trimming. Can be matrix, otuTable, or higher-
#' order phyloseq classes that contain an otuTable.
#'
#' @param flist An enclosure object, typically created with \code{\link{filterfunSample}}
#'
#' @param A An integer. The number of samples in which a taxa / species passed the filter
#' for it to be labeled TRUE in the output logical vector.
#'
#' @return A logical vector with names equal to species.names (or rownames, if matrix).
#'
#' @seealso \code{\link[genefilter]{genefilter}}, \code{\link{filterfunSample}},
#'  \code{\link[phyloseq]{t}},
#'  \code{\link{prune_species}}
#' @keywords agglomerate OTU cluster tree
#'
#' @rdname genefilterSample-methods
#' @docType methods
#' @export
#'
#' @examples #
#' ## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' ## f1  <- filterfunSample(topk(2))
#' ## wh1 <- genefilterSample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_species(wh1, testOTU)
#' ## prune_species(wh2, testOTU)
#' ## 
#' ## taxtab1 <- taxTab(matrix("abc", 5, 5))
#' ## prune_species(wh1, taxtab1)
#' ## prune_species(wh2, taxtab1)
setGeneric("genefilterSample", function(X, flist, A=1) standardGeneric("genefilterSample"))
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,matrix-method
setMethod("genefilterSample", signature("matrix"), function(X, flist, A=1){
	TFmat = apply(X, 2, flist)
	apply(TFmat, 1, function(x, A){sum(x) >= A}, A)
})
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,otuTable-method
setMethod("genefilterSample", signature("otuTable"), function(X, flist, A=1){
	if( speciesAreRows(X) ){
		genefilterSample(   as(X, "matrix"), flist, A)
	} else {
		genefilterSample( t(as(X, "matrix")), flist, A)
	}
})
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,phyloseq-method
setMethod("genefilterSample", signature("phyloseq"), function(X, flist, A=1){
	genefilterSample(otuTable(X), flist, A)
})
################################################################################
#' A sample-wise filter function builder, analogous to \code{\link[genefilter]{filterfun}}.
#'
#' See the \code{\link[genefilter]{filterfun}}, from the Bioconductor repository,
#' for a taxa-/gene-wise filter (and further examples).
#' 
#' @usage filterfunSample(...)
#'
#' @param ... A comma-separated list of functions.
#' 
#' @return An enclosure (function) that itself will return a logical vector, 
#'  according to the
#'  functions provided in the argument list, evaluated in order. The output of
#'  filterfunSample is appropriate for the `flist' argument to the 
#'  genefilterSample method.
#' 
#' @export
#' @seealso \code{\link[genefilter]{filterfun}}, \code{\link{genefilterSample}}
#' @examples
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' # f1  <- filterfunSample(topk(2))
#' # wh1 <- genefilterSample(testOTU, f1, A=2)
#' # wh2 <- c(T, T, T, F, F)
#' # prune_species(wh1, testOTU)
#' # prune_species(wh2, testOTU)
filterfunSample = function(...){
    flist <- list(...)
    if( length(flist) == 1 && is.list(flist[[1]])) { flist <- flist[[1]] }
    f = function(x){
    	# initialize fval (a logical vector)
    	fun  = flist[[1]]
    	fval = fun(x)
    	# check the remaining functions. Compare & logic, element-wise, each loop.
        for(fun in flist[-1]){
            fval = fval & fun(x)
		}
		return(fval)
	}
	class(f) <- "filterfun"
	return(f)
}
############################################################
#' Make filter fun. the most abundant \code{k} taxa
#'
#' @usage topk(k, na.rm=TRUE)
#'
#' @param k An integer, indicating how many of the most abundant taxa
#'  should be kept.
#' @param na.rm A logical. Should \code{NA}s be removed. Default is \code{TRUE}.
#'
#' @return Returns a function (enclosure) that will return TRUE
#'  for each element in the most abundant k values.
#'
#' @seealso \code{\link{topk}}, \code{\link{topf}},
#'  \code{\link{topp}}, \code{\link{rm_outlierf}}
#'
#' @export
#' 
#' @examples
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' # f1  <- filterfunSample(topk(2))
#' # wh1 <- genefilterSample(testOTU, f1, A=2)
#' # wh2 <- c(T, T, T, F, F)
#' # prune_species(wh1, testOTU)
#' # prune_species(wh2, testOTU)
topk = function(k, na.rm=TRUE){
    function(x){
		if(na.rm){x = x[!is.na(x)]}
		x >= sort(x, decreasing=TRUE)[k]
    }
}
############################################################
#' Make filter fun. that returns the most abundant \code{p} fraction of taxa
#'
#' @usage topp(p, na.rm=TRUE)
#'
#' @param p A numeric of length 1, indicating what fraction of the most abundant taxa
#'  should be kept.
#' @param na.rm A logical. Should \code{NA}s be removed. Default is \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}},
#'  that will return \code{TRUE}
#'  for each element in the most abundant p fraction of taxa.
#'
#' @seealso \code{\link{topk}}, \code{\link{topf}},
#'  \code{\link{topp}}, \code{\link{rm_outlierf}}
#'
#' @export
#'
#' @examples
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' # sampleSums(testOTU)
#' # f1  <- filterfunSample(topp(0.2))
#' # (wh1 <- genefilterSample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_species(wh1, testOTU)
#' # prune_species(wh2, testOTU)
topp <- function(p, na.rm=TRUE){
    function(x){
		if(na.rm){x = x[!is.na(x)]}
		x >= sort(x, decreasing=TRUE)[ceiling(length(x)*p)]
    }
}
################################################################################
#' Make filter fun. that returns the top f fraction of taxa in a sample.
#'
#' As opposed to \code{\link{topp}}, which gives the
#' most abundant p fraction of observed taxa (richness, instead of cumulative
#' abundance. Said another way, topf ensures a certain
#' fraction of the total sequences are retained, while topp ensures
#' that a certain fraction of taxa/species/OTUs are retained.
#'
#' @usage topf(f, na.rm=TRUE)
#' @param f Single numeric value between 0 and 1.
#' @param na.rm Logical. Should we remove NA values. Default \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}},
#'  that will return \code{TRUE}
#'  for each element in the taxa comprising the most abundant f fraction of individuals.
#'
#' @seealso \code{\link{topk}}, \code{\link{topf}},
#'  \code{\link{topp}}, \code{\link{rm_outlierf}}
#'
#' @export
#' 
#' @examples
#' # t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
#' # topf(0.6)(t1)
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' # f1  <- filterfunSample(topf(0.4))
#' # (wh1 <- genefilterSample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_species(wh1, testOTU)
#' # prune_species(wh2, testOTU)
topf <- function(f, na.rm=TRUE){
    function(x){
        if (na.rm){
            x = x[!is.na(x)]
        }
        y <- sort(x, decreasing = TRUE)
        y <- cumsum(y)/sum(x)
        return( (y <= f)[names(x)] )
    }
}
################################################################################
#' Set to FALSE any outlier species greater than f fractional abundance.
#'
#' This is for removing overly-abundant outlier taxa, not for trimming low-abundance
#' taxa.
#'
#' @usage rm_outlierf(f, na.rm=TRUE)
#'
#' @param f Single numeric value between 0 and 1. The maximum fractional abundance
#'  value that a taxa will be allowed to have in a sample without being marked
#'  for trimming.
#'
#' @param na.rm Logical. Should we remove NA values. Default \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}}.
#'
#' @seealso \code{\link{topk}}, \code{\link{topf}},
#'  \code{\link{topp}}, \code{\link{rm_outlierf}}
#'
#' @export
#' @examples
#' t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
#' rm_outlierf(0.15)(t1)
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' # speciesSums(testOTU)
#' # f1  <- filterfunSample(rm_outlierf(0.1))
#' # (wh1 <- genefilterSample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_species(wh1, testOTU)
#' # prune_species(wh2, testOTU) 
rm_outlierf <- function(f, na.rm=TRUE){
	function(x){
		if(na.rm){
			x = x[!is.na(x)]
		}
		y <- x / sum(x)
        return( y < f )
    }
}
################################################################################