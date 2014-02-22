################################################################################
# Function to create subsampled dataset 
# in which each sample has same number of total observations/counts/reads
# Note that the subsampling is random, so some noise is introduced making the
# relative abundances slightly different
################################################################################
#' Resample an OTU table such that all samples have the same library size.
#' 
#' Please note that the authors of phyloseq do not advocate using this
#' as a normalization procedure, despite its recent popularity.
#' Our justifications for using alternative approaches to address
#' disparities in library sizes have been made available as 
#' \href{http://arxiv.org/abs/1310.0424}{an article in the quantitative-biology sub-arXiv}
#' ahead of peer-reviewed publication.
#' See \code{\link{phyloseq_to_deseq2}} for a recommended alternative to rarefying.
#' Nevertheless, for comparison and demonstration, the rarefying procedure is implemented
#' here in good faith and with options we hope are useful.
#' This function uses the standard R \code{\link{sample}} function to 
#' resample from the abundance values 
#' in the \code{otu_table} component of the first argument,
#' \code{physeq}.
#' Often one of the major goals of this procedure is to achieve parity in 
#' total number of counts between samples, as an alternative to other formal
#' normalization procedures, which is why a single value for the 
#' \code{sample.size} is expected.
#' This kind of resampling can be performed with and without replacement,
#' with replacement being the more computationally-efficient, default setting.
#' See the \code{replace} parameter documentation for more details.
#' We recommended that you explicitly select a random number generator seed
#' before invoking this function, or, alternatively, that you 
#' explicitly provide a single positive integer argument as \code{rngseed}.
#'
#' This approach is sometimes mistakenly called "rarefaction", which
#' \href{http://en.wikipedia.org/wiki/Rarefaction}{in physics refers to a form of wave decompression;}
#' but in this context, ecology, the term refers to a
#' \href{http://en.wikipedia.org/wiki/Rarefaction_(ecology)}{repeated sampling procedure to assess species richness},
#' first proposed in 1968 by Howard Sanders.
#' In contrast, the procedure implemented here is used as an \emph{ad hoc} means to
#' normalize microbiome counts that have 
#' resulted from libraries of widely-differing sizes.
#' Here we have intentionally adopted an alternative
#' name, \code{rarefy}, that has also been used recently
#' to describe this process 
#' and, to our knowledge, not previously used in ecology.
#'
#' Make sure to use \code{\link{set.seed}} for exactly-reproducible results
#' of the random subsampling. 
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}} object that you
#'  want to trim/filter.
#'
#' @param sample.size (Optional). A single integer value equal to the number
#'  of reads being simulated, also known as the depth,
#'  and also equal to each value returned by \code{\link{sample_sums}}
#'  on the output. 
#'  
#' @param rngseed (Optional). A single integer value passed to 
#'  \code{\link{set.seed}}, which is used to fix a seed for reproducibly
#'  random number generation (in this case, reproducibly random subsampling).
#'  The default value is \code{711}. 
#'  If set to \code{FALSE}, then no fiddling with the RNG seed is performed,
#'  and it is up to the user to appropriately call \code{\link{set.seed}}
#'  beforehand to achieve reproducible results.
#'  
#' @param replace (Optional). Logical. Whether to sample with replacement
#'  (\code{TRUE}) or without replacement (\code{FALSE}).
#'  The default is with replacement (\code{replace=TRUE}).
#'  Two implications to consider are that
#'  (1) sampling with replacement is faster and more memory efficient
#'  as currently implemented; and
#'  (2), sampling with replacement means that there is a chance that the
#'  number of reads for a given OTU in a given sample could be larger
#'  than the original count value, as opposed to sampling without replacement
#'  where the original count value is the maximum possible.
#'  Prior to phyloseq package version number \code{1.5.20},
#'  this parameter did not exist and sampling with replacement was the only
#'  random subsampling implemented in the \code{rarefy_even_depth} function.
#'  Note that this default behavior was selected for computational efficiency,
#'  but differs from analogous functions in related packages
#'  (e.g. subsampling in QIIME).
#'  
#'  @param trimOTUs (Optional). Logical. Whether to trim OTUs
#'   from the dataset that are no longer observed in any sample
#'   (have a count of zero in every sample). 
#'   The number of OTUs trimmed, if any, is printed to
#'   standard out as a reminder.
#'   
#' @param verbose (Optional). Logical. Default is \code{TRUE}.
#'  If \code{TRUE}, extra non-warning, non-error messages are printed
#'  to standard out, describing steps in the rarefying process, 
#'  the OTUs and samples removed, etc. This can be useful the
#'  first few times the function is executed, but can be set
#'  to \code{FALSE} as-needed once behavior has been verified
#'  as expected.  
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
#' # Test with esophagus dataset
#' data("esophagus")
#' esorepT = rarefy_even_depth(esophagus, replace=TRUE)
#' esorepF = rarefy_even_depth(esophagus, replace=FALSE)
#' sample_sums(esophagus)
#' sample_sums(esorepT)
#' sample_sums(esorepF)
#' ## NRun Manually: Too slow!
#' # data("GlobalPatterns")
#' # GPrepT = rarefy_even_depth(GlobalPatterns, 1E5, replace=TRUE)
#' ## Actually just this one is slow
#' # system.time(GPrepF <- rarefy_even_depth(GlobalPatterns, 1E5, replace=FALSE))
rarefy_even_depth <- function(physeq, sample.size=min(sample_sums(physeq)),
															rngseed=FALSE, replace=TRUE, trimOTUs=TRUE, verbose=TRUE){
	
	if( as(rngseed, "logical") ){
		# Now call the set.seed using the value expected in phyloseq
		set.seed(rngseed)
    if(verbose){
      # Print to screen this value
      message("`set.seed(", rngseed, ")` was used to initialize repeatable random subsampling.")
      message("Please record this for your records so others can reproduce.")
      message("Try `set.seed(", rngseed,"); .Random.seed` for the full vector", sep="")		
      message("...")      
    }
	} else if(verbose){
		message("You set `rngseed` to FALSE. Make sure you've set & recorded\n",
				" the random seed of your session for reproducibility.\n",
				"See `?set.seed`\n")
		message("...")
	}
	
	# Make sure sample.size is of length 1.
	if( length(sample.size) > 1 ){
		warning("`sample.size` had more than one value. ", 
            "Using only the first. \n ... \n")
		sample.size <- sample.size[1]	
	}
	
	if( sample.size <= 0 ){
		stop("sample.size less than or equal to zero. ", 
         "Need positive sample size to work.")
	}
	
	# Instead of warning, expected behavior now is to prune samples
	# that have fewer reads than `sample.size`
	if( min(sample_sums(physeq)) < sample.size ){
		rmsamples = sample_names(physeq)[sample_sums(physeq) < sample.size]
    if(verbose){
      message(length(rmsamples), " samples removed",
          "because they contained fewer reads than `sample.size`.")
      message("Up to first five removed samples are: \n")
      message(rmsamples[1:min(5, length(rmsamples))], sep="\t")
      message("...")      
    }
		# Now done with notifying user of pruning, actually prune.
		physeq = prune_samples(setdiff(sample_names(physeq), rmsamples), physeq)
	}
	# initialize the subsamples phyloseq instance, newsub
	newsub <- physeq
	# enforce orientation as species-are-rows, for assignment
	if(!taxa_are_rows(newsub)){newsub <- t(newsub)}
	# apply through each sample, and replace
	newotu <- apply(otu_table(newsub), 2, rarefaction_subsample,
									sample.size=sample.size, replace=replace)
	# Add OTU names to the row indices
	rownames(newotu) <- taxa_names(physeq)
	# replace the otu_table.
	otu_table(newsub) <- otu_table(newotu, TRUE)
  if(trimOTUs){
    # Check for and remove empty OTUs
    # 1. Notify user of empty OTUs being cut.
    # 2. Cut empty OTUs
    rmtaxa = taxa_names(newsub)[taxa_sums(newsub) <= 0]
    if( length(rmtaxa) > 0 ){
      if(verbose){
        message(length(rmtaxa), "OTUs were removed because they are no longer \n",
            "present in any sample after random subsampling\n")
        message("...")
      }
      newsub = prune_taxa(setdiff(taxa_names(newsub), rmtaxa), newsub)
    }
  }
	return(newsub)
}
################################################################################
# rarefaction subsample function, one sample
################################################################################
#' @keywords internal
rarefaction_subsample <- function(x, sample.size, replace=FALSE){
	# This is a test
	# x = sample(10, 10)
	# x = 1:10
	# sample.size = 50
	#system.time(obsvec <- foreach(OTUi=1:length(x), times=x, .combine=c) %do% {rep(OTUi, times)})
	# data("GlobalPatterns")
	# sample.size = sample_sums(GlobalPatterns)[which.min(sample_sums(GlobalPatterns))]
	# x = get_taxa(GlobalPatterns, which.max(sample_sums(GlobalPatterns)))
	# Create replacement species vector
	rarvec <- numeric(length(x))	
	# Perform the sub-sampling. Suppress warnings due to old R compat issue.
	# Also, make sure to avoid errors from x summing to zero, 
  # and there are no observations to sample.
	# The initialization of rarvec above is already sufficient.
	if(sum(x) <= 0){
		# Protect against, and quickly return an empty vector, 
		# if x is already an empty count vector
		return(rarvec)
	}
	if(replace){
		# resample with replacement
		suppressWarnings(subsample <- sample(1:length(x), sample.size, replace=TRUE, prob=x))
	} else {
		# resample without replacement
		obsvec <- apply(data.frame(OTUi=1:length(x), times=x), 1, function(x){
			rep_len(x["OTUi"], x["times"])
		})
		obsvec <- unlist(obsvec, use.names=FALSE)
		# use `sample` for subsampling. Hope that obsvec doesn't overflow.
		suppressWarnings(subsample <- sample(obsvec, sample.size, replace=FALSE))
	}
	# Tabulate the results (these are already named by the order in `x`)
	sstab <- table(subsample)
	# Assign the tabulated random subsample values to the species vector
	rarvec[as(names(sstab), "integer")] <- sstab
	# Return abundance vector. Let replacement happen elsewhere.
	return(rarvec)
}
################################################################################
#' Agglomerate closely-related taxa using single-linkage clustering.
#' 
#' All tips of the tree separated by a cophenetic distance smaller than 
#' \code{h} will be agglomerated into one taxa using \code{\link{merge_taxa}}.
#' 
#' Can be used to create a non-trivial OTU Table, if a phylogenetic tree is available.
#'
#' For now, a simple, ``greedy'', single-linkage clustering is used. In future releases
#' it should be possible to specify different clustering approaches available in \code{R},
#' in particular, complete-linkage clustering appears to be used more commonly for OTU
#' clustering applications.
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}},
#'  containing a phylogenetic tree. 
#'  Alternatively, a phylogenetic tree \code{\link[ape]{phylo}} will also work.
#'
#' @param h (Optional). Numeric scalar of the height where the tree should be cut.
#' This refers to the tree resulting from hierarchical clustering
#' of \code{\link[ape]{cophenetic.phylo}(phy_tree(physeq))},
#' not necessarily the original phylogenetic tree, \code{phy_tree(physeq)}.
#' Default value is \code{0.2}. 
#' Note that this argument used to be named \code{speciationMinLength},
#' before this function/method was rewritten.
#' 
#' @param hcfun (Optional). A function. 
#'  The (agglomerative, hierarchical) clustering function to use.
#'  Good examples are
#'  \code{\link[cluster]{agnes}} and \code{\link[stats]{hclust}}.
#'  The default is \code{\link[cluster]{agnes}}.
#'  
#' @param ... (Optional). Additional named arguments to pass
#'  to \code{hcfun}. 
#'
#' @return An instance of the \code{\link{phyloseq-class}}.
#'  Or alternatively, a \code{\link{phylo}} object if the
#'  \code{physeq} argument was just a tree.
#'  In the expected-use case, the number of OTUs will be fewer
#'  (see \code{\link{ntaxa}}),
#'  after merging OTUs that are related enough to be called
#'  the same OTU. 
#'
#' @seealso 
#' 
#' \code{\link{merge_taxa}}
#' 
#' \code{\link[cluster]{agnes}}
#' 
#' \code{\link[stats]{hclust}}
#' 
#' \code{\link[ape]{cophenetic.phylo}}
#' 
#' \code{\link[ape]{phylo}}
#'
#' @importFrom cluster agnes
#'
#' @export
#'
#' @examples 
#' data("esophagus")
#' # for speed
#' esophagus = prune_taxa(taxa_names(esophagus)[1:25], esophagus)
#' plot_tree(esophagus, label.tips="taxa_names", size="abundance", title="Before tip_glom()")
#' plot_tree(tip_glom(esophagus, h=0.2), label.tips="taxa_names", size="abundance", title="After tip_glom()")
tip_glom = function(physeq, h=0.2, hcfun=agnes, ...){
  dd = as.dist(cophenetic.phylo(phy_tree(physeq)))
  psclust = cutree(as.hclust(hcfun(dd, ...)), h=h)
  cliques = levels(factor(psclust))[tapply(psclust, factor(psclust), function(x){length(x)>1})]
  # For each clique, merge taxa in it...
  for( i in cliques){
    physeq = merge_taxa(physeq, eqtaxa=names(psclust)[psclust == i])
  }
  return(physeq)
}
################################################################################
################################################################################
#' Agglomerate taxa of the same type.
#'
#' This method merges species that have the same taxonomy at a certain 
#' taxaonomic rank. 
#' Its approach is analogous to \code{\link{tip_glom}}, but uses categorical data
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
	tax <- as(access(physeq, "tax_table"), "matrix")[, 1:CN, drop=FALSE]
	tax <- apply(tax, 1, function(i){paste(i, sep=";_;", collapse=";_;")})
	
	# Remove NAs and useless from the vector/factor for looping.
	# This does not remove the taxa that have an unknown (NA)
	# taxonomic designation at this particular taxonomic rank.
	tax <- tax[ !(tax %in% bad_empty) ]
	
	# Define the OTU cliques to loop through
	spCliques <- tapply(names(tax), factor(tax), list)
	
	# Successively merge taxa in physeq.
	for( i in names(spCliques)){
		physeq <- merge_taxa(physeq, spCliques[[i]])
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
#' \href{http://cran.at.r-project.org/web/packages/picante/index.html}{prune.sample}.
#'
#' @usage prune_taxa(taxa, x)
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
#'  
#'  \code{\link{prune_samples}}
#'
#'  \href{http://cran.at.r-project.org/web/packages/picante/index.html}{prune.sample}
#'
#' @rdname prune_taxa-methods
#' @export
#' @examples
#' data("esophagus")
#' esophagus
#' plot(sort(taxa_sums(esophagus), TRUE), type="h", ylim=c(0, 50))
#' x1 = prune_taxa(taxa_sums(esophagus) > 10, esophagus) 
#' x2 = prune_taxa(names(sort(taxa_sums(esophagus), TRUE))[1:9], esophagus) 
#' identical(x1, x2)
setGeneric("prune_taxa", function(taxa, x) standardGeneric("prune_taxa"))
#' @aliases prune_taxa,NULL,ANY-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("NULL", "ANY"), function(taxa, x){
	return(x)
})
# Any prune_taxa call w/ signature starting with a logical
# converts the logical to a character vector, and then dispatches
# to more specific method.
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
#' @importFrom ape drop.tip
#' @aliases prune_taxa,character,phylo-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "phylo"), function(taxa, x){
	if( length(taxa) <= 1 ){
		# Can't have a tree with 1 or fewer tips
		warning("prune_taxa attempted to reduce tree to 1 or fewer tips.\n tree replaced with NULL.")
		return(NULL)
	} else if( setequal(taxa, taxa_names(x)) ){
		return(x)
	} else {
		return( drop.tip(x, setdiff(taxa_names(x), taxa)) )		
	}
})
#' @aliases prune_taxa,character,otu_table-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "otu_table"), function(taxa, x){
	if( setequal(taxa, taxa_names(x)) ){
		return(x)
	} else {
		taxa = intersect( taxa, taxa_names(x) )
		if( taxa_are_rows(x) ){
			return(x[taxa, , drop=FALSE])
		} else {
			return(x[, taxa, drop=FALSE])
		}
	}
})
#' @aliases prune_taxa,character,sample_data-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "sample_data"), function(taxa, x){
	return(x)
})
#' @aliases prune_taxa,character,phyloseq-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "phyloseq"), function(taxa, x){
	# Re-define `taxa` as the intersection of OTU names for each component AND `taxa`
	taxa = intersect(intersect_taxa(x), taxa)
	# Now prune them all.
	# All phyloseq objects have an otu_table slot, no need to test for existence.
	x@otu_table     = prune_taxa(taxa, otu_table(x))
	# Test if slot is present. If so, perform the component prune.
	if( !is.null(x@tax_table) ){
		x@tax_table = prune_taxa(taxa, tax_table(x))
	}
	if( !is.null(x@phy_tree) ){
		x@phy_tree  = prune_taxa(taxa, phy_tree(x))
	}
	if( !is.null(x@refseq) ){
		x@refseq    = prune_taxa(taxa, refseq(x))
	}	
	# Force index order after pruning to be the same,
	# according to the same rules as in the constructor, phyloseq()
	x = index_reorder(x, index_type="taxa")
	return(x)
})
#' @aliases prune_taxa,character,taxonomyTable-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "taxonomyTable"), function(taxa, x){
	if( setequal(taxa, taxa_names(x)) ){
		return(x)
	} else {
		taxa = intersect( taxa, taxa_names(x) )
		return( x[taxa, , drop=FALSE] )
	}
})
#' @importClassesFrom Biostrings XStringSet
#' @aliases prune_taxa,character,XStringSet-method
#' @rdname prune_taxa-methods
setMethod("prune_taxa", signature("character", "XStringSet"), function(taxa, x){
	if( setequal(taxa, taxa_names(x)) ){
		# Nothing to do, return x as-is.
		return(x)
	} else if( length(intersect(taxa, taxa_names(x))) == 0 ){
		# Informative error if intersection is zero.
		stop("prune_taxa,XStringSet: taxa and taxa_names(x) do not overlap.")		
	} else {
		# Pop the OTUs that are not in `taxa`, without reordering.
		return(x[-which(!taxa_names(x) %in% taxa)])
	}
})
################################################################################
################################################################################
#' Define a subset of samples to keep in a phyloseq object.
#' 
#' An S4 Generic method for pruning/filtering unwanted samples
#' by defining those you want to keep.
#'
#' @usage prune_samples(samples, x)
#'
#' @param samples (Required). A character vector of the samples in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept samples are TRUE, and length
#' is equal to the number of samples in object x. If \code{samples} is a named
#' logical, the samples retained is based on those names. Make sure they are
#' compatible with the \code{sample_names} of the object you are modifying (\code{x}). 
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
#' @examples
#'  data(GlobalPatterns)
#'  # Subset to just the Chlamydiae phylum.
#'  GP.chl <- subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
#'  # Remove the samples that have less than 20 total reads from Chlamydiae
#'  GP.chl <- prune_samples(sample_sums(GP.chl)>=20, GP.chl)
#'  # (p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))
setGeneric("prune_samples", function(samples, x) standardGeneric("prune_samples"))
#' @aliases prune_samples,character,otu_table-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "otu_table"), function(samples, x){
	if( setequal(samples, sample_names(x)) ){
		# If the sets of `samples` and sample_names are the same, return as-is.
		return(x)
	} else {
		samples = intersect(samples, sample_names(x))
		if( taxa_are_rows(x) ){
			return( x[, samples] )
		} else {
			return( x[samples, ] )
		}
	}
})
#' @aliases prune_samples,character,sample_data-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "sample_data"), function(samples, x){
	if( setequal(samples, sample_names(x)) ){
		# If the sets of `samples` and sample_names are the same, return as-is.
		return(x)
	} else {
		samples = intersect(samples, sample_names(x))	
		return(x[samples, ])
	}
})
#' @aliases prune_samples,character,phyloseq-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "phyloseq"), function(samples, x){
	# Re-define `samples` as the intersection of samples names for each component AND `samples`
	samples = intersect(intersect_samples(x), samples)
	# Now prune each component.
	# All phyloseq objects have an otu_table slot, no need to test for existence.
	x@otu_table = prune_samples(samples, otu_table(x))	
	if( !is.null(x@sam_data) ){
		# protect missing sample_data component. Don't need to prune if empty
		x@sam_data = prune_samples(samples, sample_data(x))
	}
	# Force sample index order after pruning to be the same,
	# according to the same rules as in the constructor, phyloseq()
	x = index_reorder(x, index_type="samples")
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
#' data(esophagus)
#' x1 = transform_sample_counts(esophagus, threshrankfun(50))
#' otu_table(x1)
#' x2 = transform_sample_counts(esophagus, rank)
#' otu_table(x2)
#' identical(x1, x2)
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
#' Transform abundance data in an \code{otu_table}, sample-by-sample.
#' 
#' This function transforms the sample counts of a taxa
#' abundance matrix according to a user-provided function.
#' The counts of each sample will be transformed individually. No sample-sample 
#' interaction/comparison is possible by this method. 
#'
#' @usage transform_sample_counts(physeq, fun, ...)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} of \code{\link{otu_table-class}}.
#'
#' @param fun (Required). A single-argument function that will be applied
#'  to the abundance counts of each sample. 
#'  Can be an anonymous \code{\link[base]{function}}.
#'  
#' @param ... (Optional). Additional, optionally-named, arguments passed to  
#'  \code{fun} during transformation of abundance data.
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
#' data(esophagus)
#' x1 = transform_sample_counts(esophagus, threshrankfun(50))
#' head(otu_table(x1), 10)
#' x2 = transform_sample_counts(esophagus, rank)
#' head(otu_table(x2), 10)
#' identical(x1, x2)
#' x3 = otu_table(esophagus) + 5
#' x3 = transform_sample_counts(x3, log)
#' head(otu_table(x3), 10)
#' x4 = transform_sample_counts(esophagus, function(x) round(x^2.2, 0))
#' head(otu_table(x4), 10)
transform_sample_counts <- function(physeq, fun, ...){
	# Test the user-provided function returns a vector of the same length as input.
	if( !identical(length(fun(1:10)), 10L) ){stop("`fun` not valid function.")}
	# Check orientation, transpose if-needed to make apply work properly.
	if( taxa_are_rows(physeq) ){
		newphyseq = apply(as(otu_table(physeq), "matrix"), 2, fun, ...)
		if( identical(ntaxa(physeq), 1L) ){
			# Fix the dropped index when only 1 OTU.
			newphyseq <- matrix(newphyseq, 1L, nsamples(physeq), TRUE,
													list(taxa_names(physeq), sample_names(physeq)))
		}
	} else {
		newphyseq = apply(t(as(otu_table(physeq), "matrix")), 2, fun, ...)
		if( identical(ntaxa(physeq), 1L) ){
			# Fix the dropped index when only 1 OTU.
			newphyseq <- matrix(newphyseq, 1L, nsamples(physeq), TRUE,
													list(taxa_names(physeq), sample_names(physeq)))
		}
		newphyseq = t(newphyseq)
	}
	# Check that original and new dimensions agree. Error if not.
	if( !identical(dim(newphyseq), dim(otu_table(physeq))) ){
		stop("Dimensions of OTU table change after apply-ing function. \n",
			"       Please check both function and table")
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
#' ## wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
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
#' A sample-wise filter function builder
#' analogous to \code{\link[genefilter]{filterfun}}.
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
#' # Use simulated abundance matrix
#' set.seed(711)
#' testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' f1  <- filterfun_sample(topk(2))
#' wh1 <- genefilter_sample(testOTU, f1, A=2)
#' wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' prune_taxa(wh1, testOTU)
#' prune_taxa(wh2, testOTU)
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
#' set.seed(711)
#' testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' f1  <- filterfun_sample(topk(2))
#' wh1 <- genefilter_sample(testOTU, f1, A=2)
#' wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' prune_taxa(wh1, testOTU)
#' prune_taxa(wh2, testOTU)
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
#' set.seed(711)
#' testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' sample_sums(testOTU)
#' f1  <- filterfun_sample(topp(0.2))
#' (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' prune_taxa(wh1, testOTU)
#' prune_taxa(wh2, testOTU)
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
#' t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
#' topf(0.6)(t1)
#' ## Use simulated abundance matrix
#' set.seed(711)
#' testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' f1  <- filterfun_sample(topf(0.4))
#' (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' prune_taxa(wh1, testOTU)
#' prune_taxa(wh2, testOTU)
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
#' set.seed(711)
#' testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#' taxa_sums(testOTU)
#' f1  <- filterfun_sample(rm_outlierf(0.1))
#' (wh1 <- genefilter_sample(testOTU, f1, A=1))
#' wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' prune_taxa(wh1, testOTU)
#' prune_taxa(wh2, testOTU) 
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
