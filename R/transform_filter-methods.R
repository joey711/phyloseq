################################################################################
# Function to create subsampled dataset 
# in which each sample has same number of total observations/counts/reads
# Note that the subsampling is random, so some noise is introduced making the
# relative abundances slightly different
################################################################################
#' Perform a random subsampling of an OTU table to a level of even depth.
#' 
#' This function uses the \code{\link{sample}} function to randomly subset from the 
#' abundance values in each sample of the \code{otu_table} component in the
#' \code{physeq} argument.
#' Sampling is performed with replacement from a vector of taxa indices,
#' with length equal to the argument to \code{sample.size},
#' and probability according to the abundances for that sample in \code{physeq}.
#'
#' This is sometimes (somewhat mistakenly) called "rarefaction", 
#' though it actually a single random subsampling procedure in this case. 
#' The original rarefaction procedure includes many
#' random subsampling iterations at increasing depth as a means to
#' infer richness/alpha-diversity
#'
#' Make sure to use \code{\link{set.seed}} for exactly-reproducible results
#' of the random subsampling. 
#'
#' @usage rarefy_even_depth(physeq, sample.size=min(sample_sums(physeq)))
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}} object that you
#'  want to trim/filter.
#'
#' @param sample.size (Optional). A single integer value equal to the number
#'  of reads being simulated, also known as the depth,
#'  and also equal to each value returned by \code{\link{sample_sums}}
#'  on the output. 
#'
#' @return An object of class \code{phyloseq}. 
#' Only the \code{otu_table} component is modified.
#'
#' @seealso
#' \code{\link{sample}}
#' 
#' \code{\link{set.seed}}
#'
#' @export
#'
#' @examples
#' set.seed(711)
#' # Test with esophagus dataset
#' data("esophagus")
#' eso <- rarefy_even_depth(esophagus)
#' plot(as(otu_table(eso), "vector"), as(otu_table(esophagus), "vector"))
#' UniFrac(eso); UniFrac(esophagus)
#' # Test with GlobalPatterns dataset
#' data("GlobalPatterns")
#' GP.chl <- subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
#' # remove the samples that have less than 20 total reads from Chlamydiae
#' GP.chl <- prune_samples(names(which(sample_sums(GP.chl)>=20)), GP.chl)
#' # # (p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))
#' GP.chl.r <- rarefy_even_depth(GP.chl)
#' plot(as(otu_table(GP.chl.r), "vector"), as(otu_table(GP.chl), "vector"))
#' # Try ordination of GP.chl and GP.chl.r (default distance is unweighted UniFrac)
#' plot_ordination(GP.chl, ordinate(GP.chl, "MDS"), color="SampleType") #+ geom_point(size=5)
#' plot_ordination(GP.chl.r, ordinate(GP.chl.r, "MDS"), color="SampleType") #+ geom_point(size=5)
rarefy_even_depth <- function(physeq, sample.size=min(sample_sums(physeq))){
	sample.size <- sample.size[1]
	
	if( sample.size <= 0 ){
		stop("sample.size less than or equal to zero. Need positive sample size to work.")
	}
	if( min(sample_sums(physeq)) < sample.size ){
		warning("Strange behavior expected for samples with fewer observations than sample.size")
	}
	# initialize the subsamples phyloseq instance, newsub
	newsub <- physeq
	# enforce orientation as species-are-rows, for assignment
	if(!taxa_are_rows(newsub)){newsub <- t(newsub)}
	# apply through each sample, and replace
	newotu <- apply(otu_table(newsub), 2, rarefaction_subsample, sample.size)
	# Add species names to the row indices
	rownames(newotu) <- taxa_names(physeq)
	# replace the otu_table.
	otu_table(newsub) <- otu_table(newotu, TRUE)
	return(newsub)
}
################################################################################
# rarefaction subsample function, one sample
################################################################################
#' @keywords internal
rarefaction_subsample <- function(x, sample.size){
	# Create replacement species vector
	rarvec <- numeric(length(x))	
	# Perform the subsampling. Suppress warnings due to old R compat issue.
	# Also, make sure to avoid errors from x summing to zero, and there are no observations to sample.
	# The initialization of rarvec above is already sufficient.
	if(sum(x) > 0){
		suppressWarnings(subsample <- sample(1:length(x), sample.size, TRUE, prob=x))
		# Tabulate the results
		sstab <- table(subsample)
		# Assign the tabulated random subsample values to the species vector
		rarvec[as(names(sstab), "integer")] <- sstab
	}
	# Return abundance vector. Let replacement happen elsewhere.
	return(rarvec)
}
################################################################################
#' Agglomerate closely-related taxa using single-linkage clustering.
#' 
#' All tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated into one taxa using \code{merge_taxa}.
#' 
#' Can be used to create a non-trivial OTU Table, if a phylogenetic tree is available.
#'
#' For now, a simple, ``greedy'', single-linkage clustering is used. In future releases
#' it should be possible to specify different clustering approaches available in \code{R},
#' in particular, complete-linkage clustering appears to be used more commonly for OTU
#' clustering applications.
#'
#' @usage tip_glom(tree, OTU, speciationMinLength=0.02)
#'
#' @param tree \code{\link{phyloseq-class}}, containing an OTU Table and
#'  phylogenetic tree. If, alternatively, \code{tree} is a \code{\link{phylo-class}},
#'  then \code{OTU} is required.
#'
#' @param OTU An otu_table object. Optional. Ignored if \code{tree} is a 
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
#' which case \code{tip_glom} returns an \code{phyloseq} object.
#'
#' @rdname tip_glom-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # # data(phylocom)
#' # # # otu  <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
#' # # # x1   <- phyloseq(otu, phylocom$phylo)
#' # # # print(x1); par(mfrow=c(2, 1)); plot(phy_tree(x1))
#' # # # x2 <- tip_glom(x1, speciationMinLength = 2.5)
#' # # # plot(phy_tree(x2))
#' # # # ## Try on example datset 1
#' # # # data(GlobalPatterns); ntaxa(GlobalPatterns)
#' # # # ex7 <- tip_glom(GlobalPatterns, speciationMinLength = 0.05)
#' # # # ntaxa(ex7)
#' # data(esophagus); ntaxa(esophagus); par(mfrow=c(2, 1)); plot(phy_tree(esophagus))
#' # phy_tree(esophagus)$edge.length
#' # x3 <- tip_glom(esophagus, speciationMinLength = 0.20)
#' # ntaxa(x3); plot(phy_tree(x3))
setGeneric("tip_glom", function(tree, OTU, speciationMinLength=0.02) standardGeneric("tip_glom"))
#' @rdname tip_glom-methods
#' @aliases tip_glom,phylo,otu_table-method
setMethod("tip_glom", signature("phylo", "otu_table"), function(tree, OTU, speciationMinLength=0.02){
	# dispatch as the combined-object (phyloseq-class), auto coherence.
	tip_glom( phyloseq(OTU, tree), speciationMinLength=speciationMinLength)
})
#' @rdname tip_glom-methods
#' @aliases tip_glom,phyloseq,ANY-method
setMethod("tip_glom", signature("phyloseq"), function(tree, speciationMinLength=0.02){
	tip_glom.internal(tree, speciationMinLength=speciationMinLength)
})
#' @rdname tip_glom-methods
#' @aliases tip_glom,phylo,ANY-method
setMethod("tip_glom", signature("phylo"), function(tree, speciationMinLength=0.02){
	tip_glom.internal(tree, speciationMinLength=speciationMinLength)
})
################################################################################
#' Internal function for tiplgom.
#' 
#' Internal function, users should use the S4 method \code{\link{tip_glom}}.
#' Tree can be a \code{\link{phyloseq-class}} that contains a phylogenetic tree, 
#' This is because \code{\link{merge_taxa}} can
#' handle all the relevant objects, as can \code{\link{getTipDistMatrix}}.
#' Create Non-trivial OTU table, by agglomerating nearby tips.
#' tip_glom.internal is called by the S4 \code{tip_glom} methods. It is useful if 
#' a motivated user wants to see the internals of the implementation. By design
#' it lacks explicit object handling. Use \code{\link{tip_glom}} instead.
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
#' @seealso tip_glom
#' @importFrom igraph0 graph.adjacency
#' @importFrom igraph0 get.edgelist
#' @keywords internal
tip_glom.internal <- function(tree, speciationMinLength){
	# Create adjacency matrix, where tips are adjacent
	# if their distance is below the threshold speciationMinLength
	tipAdjacent <- (getTipDistMatrix( tree ) < speciationMinLength)
	# Define igraph0 object based on the tree-tip adjacenecy matrix
	ig          <- graph.adjacency(tipAdjacent, diag=FALSE)
	# Define the species cliques to loop through
	spCliques   <- edgelist2clique( get.edgelist(ig) )
	# successively merge
	for( i in 1:length(spCliques)){
		tree <- merge_taxa(tree, eqspecies=spCliques[[i]])
	}
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
#' @seealso tip_glom
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
	getTipDistMatrix( phy_tree(tree) )
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
#' This method merges species that have the same taxonomy at a certain 
#' taxaonomic rank. 
#' Its approach is analogous to \code{tip_glom}, but uses categorical data
#' instead of a tree. In principal, other categorical data known for all taxa
#' could also be used in place of taxonomy,
#' but for the moment, this must be stored in the \code{taxonomyTable}
#' of the data. Also, columns/ranks to the right of the rank chosen to use
#' for agglomeration will be replaced with \code{NA},
#' because they should be meaningless following agglomeration.
#'
#' @usage tax_glom(physeq, taxrank=rank_names(physeq)[1], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} or \code{\link{otu_table}}.
#'
#' @param taxrank A character string specifying the taxonomic level
#'  that you want to agglomerate over.
#'  Should be among the results of \code{rank_names(physeq)}.
#'  The default value is \code{rank_names(physeq)[1]},
#'  which may agglomerate too broadly for a given experiment.
#'  You are strongly encouraged to try different values for this argument.
#'
#' @param NArm (Optional). Logical, length equal to one. Default is \code{TRUE}.
#'  CAUTION. The decision to prune (or not) taxa for which you lack categorical
#'  data could have a large effect on downstream analysis. You may want to
#'  re-compute your analysis under both conditions, or at least think carefully
#'  about what the effect might be and the reasons explaining the absence of 
#'  information for certain taxa. In the case of taxonomy, it is often a result 
#'  of imprecision in taxonomic designation based on short phylogenetic sequences
#'  and a patchy system of nomenclature. If this seems to be an issue for your
#'  analysis, think about also trying the nomenclature-agnostic \code{\link{tip_glom}}
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
#' @seealso
#' \code{\link{tip_glom}}
#' 
#' \code{\link{prune_taxa}}
#' 
#' \code{\link{merge_taxa}}
#' 
#' @export
#'
#' @examples
#' # data(GlobalPatterns)
#' # ## print the available taxonomic ranks
#' # colnames(tax_table(GlobalPatterns))
#' # ## agglomerate at the Family taxonomic rank
#' # (x1 <- tax_glom(GlobalPatterns, taxrank="Family") )
#' # ## How many taxa before/after agglomeration?
#' # ntaxa(GlobalPatterns); ntaxa(x1)
#' # ## Look at enterotype dataset...
#' # data(enterotype)
#' # ## print the available taxonomic ranks. Shows only 1 rank available, not useful for tax_glom
#' # colnames(tax_table(enterotype))
tax_glom <- function(physeq, taxrank=rank_names(physeq)[1],
					NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){

	# Error if tax_table slot is empty
	if( is.null(access(physeq, "tax_table")) ){
		stop("The tax_glom() function requires that physeq contain a taxonomyTable")
	}
	
	# Error if bad taxrank
	if( !taxrank[1] %in% rank_names(physeq) ){
		stop("Bad taxrank argument. Must be among the values of rank_names(physeq)")		
	}

	# Make a vector from the taxonomic data.
	CN  <- which( rank_names(physeq) %in% taxrank[1] )
	tax <- as(access(physeq, "tax_table"), "matrix")[, CN]
	
	# if NArm is TRUE, remove the empty, white-space, NA values from 
	if( NArm ){
		keep_species <- names(tax)[ !(tax %in% bad_empty) ]
		physeq <- prune_taxa(keep_species, physeq)
	}

	# Concatenate data up to the taxrank column, use this for agglomeration
	tax <- as(access(physeq, "tax_table"), "matrix")[, 1:CN]
	tax <- apply(tax, 1, function(i){paste(i, sep=";_;", collapse=";_;")})
	
	# Remove NAs and useless from the vector/factor for looping.
	# This does not remove the taxa that have an unknown (NA)
	# taxonomic designation at this particular taxonomic rank.
	tax <- tax[ !(tax %in% bad_empty) ]
	
	# Define the species cliques to loop through
	spCliques <- tapply(names(tax), factor(tax), list)
	
	# Successively merge taxa in physeq.
	for( i in names(spCliques)){
		physeq <- merge_taxa(physeq, eqspecies=spCliques[[i]])
	}
	
	# "Empty" the values to the right of the rank, using NA_character_.
	if( CN < length(rank_names(physeq)) ){
		badcolumns <- (CN+1):length(rank_names(physeq))
		tax_table(physeq)[, badcolumns] <- NA_character_
	}
	
	# Return.
	return(physeq)
}
################################################################################
################################################################################
#' Prune unwanted OTUs / taxa from a phylogenetic object.
#' 
#' An S4 Generic method for removing (pruning) unwanted OTUs/taxa from phylogenetic
#' objects, including phylo-class trees, as well as native phyloseq package
#' objects. This is particularly useful for pruning a phyloseq object that has
#' more than one component that describes OTUs.
#' Credit: the \code{phylo}-class version is adapted from
#' \code{\link[picante]{prune.sample}}.
#'
#' @param taxa (Required). A character vector of the taxa in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept taxa are TRUE, and length
#' is equal to the number of taxa in object x. If \code{taxa} is a named
#' logical, the taxa retained are based on those names. Make sure they are
#' compatible with the \code{taxa_names} of the object you are modifying (\code{x}). 
#'
#' @param x (Required). A phylogenetic object, including \code{phylo} trees,
#' as well as all phyloseq classes that represent taxa. If the function
#' \code{\link{taxa_names}} returns a non-\code{NULL} value, then your object
#' can be pruned by this function.
#'
#' @return The class of the object returned by \code{prune_taxa} matches
#' the class of the argument, \code{x}.
#'
#' @seealso
#'  \code{\link{prune_taxa}}
#'
#'  \code{\link[picante]{prune.sample}}
#'
#' @rdname prune_taxa-methods
#' @export
#' @examples #
#' ## testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' ## f1  <- filterfun_sample(topk(2))
#' ## wh1 <- genefilter_sample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_taxa(wh1, testOTU)
#' ## prune_taxa(wh2, testOTU)
#' ## 
#' ## tax_table1 <- tax_table(matrix("abc", 5, 5))
#' ## prune_taxa(wh1, tax_table1)
#' ## prune_taxa(wh2, tax_table1)
setGeneric("prune_taxa", function(taxa, x) standardGeneric("prune_taxa"))
################################################################################
#' @aliases prune_taxa,NULL,ANY-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("NULL"), function(taxa, x){
	return(x)
})
################################################################################
# import covering ape::drop.tip
#' @import ape
#' @aliases prune_taxa,character,phylo-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "phylo"), function(taxa, x){
	if( length(taxa) <= 1 ){
		# Can't have a tree with 1 or fewer tips
		warning("prune_taxa attempted to reduce tree to 1 or fewer tips.\n tree replaced with NULL.")
		return(NULL)
	}
	trimTaxa <- setdiff(x$tip.label, taxa)
	if( length(trimTaxa) > 0 ){
		return( drop.tip(x, trimTaxa) )
	} else {
		return(x)
	}
})
################################################################################
#' @aliases prune_taxa,character,otu_table-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "otu_table"), function(taxa, x){
	taxa <- intersect( taxa, taxa_names(x) )
	if( taxa_are_rows(x) ){
		x[taxa, , drop=FALSE]
	} else {
		x[, taxa, drop=FALSE]
	}	
})
################################################################################
#' @aliases prune_taxa,character,sample_data-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "sample_data"), function(taxa, x){
	return(x)
})
################################################################################
#' @aliases prune_taxa,character,phyloseq-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "phyloseq"), 
		function(taxa, x){
			
	# Save time and return if the union of all component taxa names
	# captured by taxa_names(x) is same as taxa. 
	if( setequal(taxa_names(x), taxa) ){
		return(x)
	} else {	
		# All phyloseq objects have an otu_table slot, no need to test.
		x@otu_table   <- prune_taxa(taxa, otu_table(x))
		
		# Test if slot is present. If so, perform the component prune.
		if( !is.null(access(x, "tax_table")) ){
			x@tax_table <- prune_taxa(taxa, tax_table(x))
		}
		if( !is.null(access(x, "phy_tree")) ){
			x@phy_tree    <- prune_taxa(taxa, phy_tree(x))
		}
		return(x)
	}
})
################################################################################
#' @aliases prune_taxa,character,taxonomyTable-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "taxonomyTable"), 
		function(taxa, x){
	taxa <- intersect( taxa, taxa_names(x) )
	return( x[taxa, , drop=FALSE] )
})
################################################################################
#' @aliases prune_taxa,logical,ANY-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("logical", "ANY"), function(taxa, x){
	# Check that logical has same length as ntaxa, stop if not.
	if( !identical(length(taxa), ntaxa(x)) ){
		stop("logical argument to taxa is wrong length. Should equal ntaxa(x)")
	} else {
		# Pass on to names-based prune_taxa method
		return( prune_taxa(taxa_names(x)[taxa], x) )		
	}
})
################################################################################
################################################################################
#' Prune unwanted samples from a phyloseq object.
#' 
#' An S4 Generic method for removing (pruning) unwanted samples.
#'
#' @usage prune_samples(samples, x)
#'
#' @param samples (Required). A character vector of the samples in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept samples are TRUE, and length
#' is equal to the number of samples in object x. If \code{samples} is a named
#' logical, the samples retained is based on those names. Make sure they are
#' compatible with the \code{taxa_names} of the object you are modifying (\code{x}). 
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
#'  data(GlobalPatterns)
#'  # Subset to just the Chlamydiae phylum.
#'  GP.chl <- subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
#'  # Remove the samples that have less than 20 total reads from Chlamydiae
#'  GP.chl <- prune_samples(sampleSums(GP.chl)>=20, GP.chl)
#'  # (p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))
setGeneric("prune_samples", function(samples, x) standardGeneric("prune_samples"))
#' @aliases prune_samples,character,otu_table-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "otu_table"), function(samples, x){
	if( taxa_are_rows(x) ){
		x[, samples]
	} else {
		x[samples, ]
	}
})
#' @aliases prune_samples,character,sample_data-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "sample_data"), function(samples, x){
	x[samples, ]
})
#' @aliases prune_samples,character,phyloseq-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "phyloseq"), function(samples, x){
	# protect missing sample_data component. Don't need to prune if empty
	if( !is.null(access(x, "sam_data", FALSE)) ){
		x@sam_data  <- prune_samples(samples, access(x, "sam_data", FALSE) )
	}
	# Don't need to protect otu_table, it is mandatory for phyloseq-class
	x@otu_table <- prune_samples(samples, access(x, "otu_table", FALSE) )
	return(x)
})
# A logical should specify the samples to keep, or not. Have same length as nsamples(x) 
#' @aliases prune_samples,logical,ANY-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("logical", "ANY"), function(samples, x){
	# Check that logical has same length as nsamples, stop if not.
	if( !identical(length(samples), nsamples(x)) ){
		stop("logical argument to samples is wrong length. Should equal nsamples(x)")
	} else {
		# Pass on to names-based prune_samples method
		return( prune_samples(sample_names(x)[samples], x) )		
	}
})
################################################################################
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
#' @seealso \code{\link{transform_sample_counts}}, \code{\link{rank}}, \code{\link{threshrankfun}}
#' @export 
#' @examples #
#' (a_vector <- sample(0:10, 100, TRUE))
#' threshrank(a_vector, 5, keep0s=TRUE)
#' data(GlobalPatterns)
#' GP <- GlobalPatterns
#' ## These three approaches result in identical otu_table
#' (x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
#' (x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
#' identical(x1, x2)
#' (x3 <- otu_table(apply(otu_table(GP), 2, threshrank, thresh=500), taxa_are_rows(GP)) )
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
#' like \code{\link{transform_sample_counts}}.
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
#' @seealso \code{\link{transform_sample_counts}}, \code{\link{threshrankfun}},
#'  \code{\link{threshrank}}
#' @export
#' @examples
#' data(GlobalPatterns)
#' GP <- GlobalPatterns
#' ## These three approaches result in identical otu_table
#' (x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
#' (x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
#' identical(x1, x2)
#' (x3 <- otu_table(apply(otu_table(GP), 2, threshrank, thresh=500), taxa_are_rows(GP)) )
#' identical(x1, x3)
threshrankfun <- function(thresh, keep0s=FALSE, ...){
	function(x){
		threshrank(x, thresh, keep0s=FALSE, ...)
	}
}
################################################################################
#' Transpose \code{\link{otu_table-class}} or \code{\link{phyloseq-class}}
#'
#' Extends the base transpose method, \code{\link[base]{t}}.
#'
#' @usage t(x)
#'
#' @param x An \code{otu_table} or \code{\link{phyloseq-class}}.
#'
#' @return The class of the object returned by \code{t} matches
#' the class of the argument, \code{x}. The \code{otu_table} is
#' transposed, and \code{\link{taxa_are_rows}} value is toggled.
#'
#' @name t
#' @rdname transpose-methods
#' @docType methods
#' @export
#' @examples
#' data(GlobalPatterns)
#' otu_table(GlobalPatterns)
#' t( otu_table(GlobalPatterns) )
setGeneric("t")
#' @aliases t,otu_table-method
#' @rdname transpose-methods
setMethod("t", signature("otu_table"), function(x){
	#new("otu_table", t(x@.Data), taxa_are_rows = (!taxa_are_rows(x)))
	x <- otu_table( t(as(x, "matrix")), taxa_are_rows=(!taxa_are_rows(x)) )
	return(x)
})
################################################################################
#' @aliases t,phyloseq-method
#' @rdname transpose-methods
setMethod("t", signature("phyloseq"), function(x){
	x@otu_table <- t( otu_table(x) )
	return(x)
})
################################################################################
#' Transform the abundance count data in an \code{otu_table}, sample-by-sample.
#' 
#' This function transforms the sample counts of a taxa
#' abundance matrix according to a user-provided function.
#' The counts of each sample will be transformed individually. No sample-sample 
#' interaction/comparison is possible by this method. 
#'
#' @usage transform_sample_counts(physeq, fun)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} of \code{\link{otu_table-class}}.
#'
#' @param fun (Required). A single-argument function that will be applied
#'  to the abundance counts of each sample. Can be an anonymous \code{\link[base]{function}}.
#' 
#' @return A transformed \code{otu_table} -- or \code{phyloseq} object with its
#'  transformed \code{otu_table}. 
#'  In general, trimming is not expected by this 
#'  method, so it is suggested that the user provide only functions that return
#'  a full-length vector. Filtering/trimming can follow, for which the 
#'  \code{\link{genefilter_sample}} and \code{\link{prune_taxa}} functions
#'  are suggested.
#'
#' @seealso \code{\link{threshrankfun}}, \code{\link{rank}}, \code{\link{log}}
#'
#' @docType methods
#' @aliases transform_sample_counts transformSampleCounts
#' @rdname transformcounts
#' @export
#'
#' @examples #
#' data(GlobalPatterns)
#' GP <- GlobalPatterns
#' ## transform_sample_counts can work on phyloseq-class, modifying otu_table only
#' (GPr <- transform_sample_counts(GP, rank) )
#' ## These two approaches result in identical otu_table
#' (x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
#' (x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
#' identical(x1, x2)
transform_sample_counts <- function(physeq, fun){
	if( taxa_are_rows(physeq) ){
		newphyseq <- apply(as(otu_table(physeq), "matrix"), 2, fun)
	} else {
		newphyseq <- apply(as(otu_table(physeq), "matrix"), 1, fun)
	}
	otu_table(physeq) <- otu_table(newphyseq, taxa_are_rows=taxa_are_rows(physeq))
	return(physeq)
}
####################################################################################
#' @rdname transformcounts
#' @export
transformSampleCounts <- transform_sample_counts
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
#' logical having length equal to nrow()/\code{\link{ntaxa}} is returned, indicating which
#' should be kept. This output can be provided
#' directly to OTU trimming function, \code{\link{prune_taxa}}.
#' By contrast, \code{\link[genefilter]{genefilter}}, 
#' of the genefilter package in Bioconductor,
#' works only on the rows of a matrix. Note that, because \code{\link{otu_table-class}}
#' inherits directly from the \code{\link{matrix-class}}, an unmodified
#' otu_table can be provided to \code{genefilter}, but be mindful of the orientation
#' of the otu_table (use \code{\link{taxa_are_rows}}),
#' and transpose (\code{\link[phyloseq]{t}}) if needed.
#'
#' @usage genefilter_sample(X, flist, A=1)
#'
#' @param X The object that needs trimming. Can be matrix, otu_table, or higher-
#' order phyloseq classes that contain an otu_table.
#'
#' @param flist An enclosure object, typically created with \code{\link{filterfun_sample}}
#'
#' @param A An integer. The number of samples in which a taxa / OTUs passed the filter
#' for it to be labeled TRUE in the output logical vector.
#'
#' @return A logical vector with names equal to taxa_names (or rownames, if matrix).
#'
#' @seealso \code{\link[genefilter]{genefilter}}, \code{\link{filterfun_sample}},
#'  \code{\link[phyloseq]{t}},
#'  \code{\link{prune_taxa}}
#' @keywords agglomerate OTU cluster tree
#'
#' @rdname genefilter_sample-methods
#' @docType methods
#' @export
#'
#' @examples #
#' ## testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' ## f1  <- filterfun_sample(topk(2))
#' ## wh1 <- genefilter_sample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_taxa(wh1, testOTU)
#' ## prune_taxa(wh2, testOTU)
#' ## 
#' ## tax_table1 <- tax_table(matrix("abc", 5, 5))
#' ## prune_taxa(wh1, tax_table1)
#' ## prune_taxa(wh2, tax_table1)
setGeneric("genefilter_sample", function(X, flist, A=1) standardGeneric("genefilter_sample"))
#' @rdname genefilter_sample-methods
#' @aliases genefilter_sample,matrix-method
setMethod("genefilter_sample", signature("matrix"), function(X, flist, A=1){
	TFmat = apply(X, 2, flist)
	apply(TFmat, 1, function(x, A){sum(x) >= A}, A)
})
#' @rdname genefilter_sample-methods
#' @aliases genefilter_sample,otu_table-method
setMethod("genefilter_sample", signature("otu_table"), function(X, flist, A=1){
	if( taxa_are_rows(X) ){
		genefilter_sample(   as(X, "matrix"), flist, A)
	} else {
		genefilter_sample( t(as(X, "matrix")), flist, A)
	}
})
#' @rdname genefilter_sample-methods
#' @aliases genefilter_sample,phyloseq-method
setMethod("genefilter_sample", signature("phyloseq"), function(X, flist, A=1){
	genefilter_sample(otu_table(X), flist, A)
})
################################################################################
#' A sample-wise filter function builder, analogous to \code{\link[genefilter]{filterfun}}.
#'
#' See the \code{\link[genefilter]{filterfun}}, from the Bioconductor repository,
#' for a taxa-/gene-wise filter (and further examples).
#' 
#' @usage filterfun_sample(...)
#'
#' @param ... A comma-separated list of functions.
#' 
#' @return An enclosure (function) that itself will return a logical vector, 
#'  according to the
#'  functions provided in the argument list, evaluated in order. The output of
#'  filterfun_sample is appropriate for the `flist' argument to the 
#'  genefilter_sample method.
#' 
#' @export
#' @seealso \code{\link[genefilter]{filterfun}}, \code{\link{genefilter_sample}}
#' @examples
#' ## Use simulated abundance matrix
#' # set.seed(711)
#' # testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' # f1  <- filterfun_sample(topk(2))
#' # wh1 <- genefilter_sample(testOTU, f1, A=2)
#' # wh2 <- c(T, T, T, F, F)
#' # prune_taxa(wh1, testOTU)
#' # prune_taxa(wh2, testOTU)
filterfun_sample = function(...){
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
################################################################################
#' Filter taxa based on across-sample OTU abundance criteria
#'
#' This function is directly analogous to the
#' \code{\link[genefilter]{genefilter}} function for microarray filtering,
#' but is used for filtering OTUs from phyloseq objects.
#' It applies an arbitrary set of functions ---
#' as a function list, for instance, created by \code{\link[genefilter]{filterfun}} ---
#' as across-sample criteria, one OTU at a time.
#' It takes as input a phyloseq object,
#' and returns a logical vector
#' indicating whether or not each OTU passed the criteria.
#' Alternatively, if the \code{"prune"} option is set to \code{FALSE},
#' it returns the already-trimmed version of the phyloseq object.
#' 
#' @usage filter_taxa(physeq, flist, prune=FALSE)
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}} object that you
#'  want to trim/filter.
#'
#' @param flist (Required). A function or list of functions that take a vector
#'  of abundance values and return a logical. Some canned useful function types
#'  are included in the \code{genefilter}-package.
#'
#' @param prune (Optional). A logical. Default \code{FALSE}. If \code{TRUE}, then
#'  the function returns the pruned \code{\link{phyloseq-class}} object, rather
#'  than the logical vector of taxa that passed the filter.
#' 
#' @return A logical vector equal to the number of taxa in \code{physeq}.
#'  This can be provided directly to \code{\link{prune_taxa}} as first argument.
#'  Alternatively, if \code{prune==TRUE}, the pruned \code{\link{phyloseq-class}} 
#'  object is returned instead.
#' 
#' @export
#' @seealso 
#' \code{\link[genefilter]{filterfun}},
#' \code{\link{genefilter_sample}},
#' \code{\link{filterfun_sample}}
#' 
#' @examples
#'  data("enterotype")
#'  require("genefilter")
#'  flist    <- filterfun(kOverA(5, 2e-05))
#'  ent.logi <- filter_taxa(enterotype, flist)
#'  ent.trim <- filter_taxa(enterotype, flist, TRUE)
#'  identical(ent.trim, prune_taxa(ent.logi, enterotype)) 
#'  identical(sum(ent.logi), ntaxa(ent.trim))
#'  filter_taxa(enterotype, flist, TRUE)
filter_taxa <- function(physeq, flist, prune=FALSE){
	# access OTU table
	OTU <- access(physeq, "otu_table", TRUE)
	# Enforce orientation (we are filtering taxa, not samples)
	if(!taxa_are_rows(OTU)) {
		OTU <- t(OTU)
	}
	# Coerce to vanilla matrix
	OTU <- as(OTU, "matrix")
	# Apply filtering function(s), get logical of length ntaxa(physeq)
	ans <- apply(OTU, 1, flist)
	# sanity check
	if( ntaxa(physeq) != length(ans) ){
		stop("Logic error in applying function(s). Logical result not same length as ntaxa(physeq)")
	}
	# Now return logical or pruned phyloseq-class instance.	
	if( prune ){
		return( prune_taxa(ans, physeq) )
	} else {
		return( ans )		
	}
}
################################################################################
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
#' # testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' # f1  <- filterfun_sample(topk(2))
#' # wh1 <- genefilter_sample(testOTU, f1, A=2)
#' # wh2 <- c(T, T, T, F, F)
#' # prune_taxa(wh1, testOTU)
#' # prune_taxa(wh2, testOTU)
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
#' @return A function (enclosure), suitable for \code{\link{filterfun_sample}},
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
#' # testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' # sample_sums(testOTU)
#' # f1  <- filterfun_sample(topp(0.2))
#' # (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_taxa(wh1, testOTU)
#' # prune_taxa(wh2, testOTU)
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
#' @return A function (enclosure), suitable for \code{\link{filterfun_sample}},
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
#' # testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' # f1  <- filterfun_sample(topf(0.4))
#' # (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_taxa(wh1, testOTU)
#' # prune_taxa(wh2, testOTU)
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
#' @return A function (enclosure), suitable for \code{\link{filterfun_sample}}.
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
#' # testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' # taxa_sums(testOTU)
#' # f1  <- filterfun_sample(rm_outlierf(0.1))
#' # (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' # wh2 <- c(T, T, T, F, F)
#' # prune_taxa(wh1, testOTU)
#' # prune_taxa(wh2, testOTU) 
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
