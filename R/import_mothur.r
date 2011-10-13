################################################################################
#' Show cutoff values available in a mothur list file
#'
#' This is a helper function to report back to the user the different cutoff
#' values available in a given \emph{list} file created by the OTU clustering
#' and analysis package called \emph{mothur}
#'
#' @usage show_mothur_list_cutoffs(mothur_list_file)
#'
#' @param mothur_list_file The list file name and/or location as produced by \emph{mothur}.
#'
#' @return A character vector of the different cutoff values contained in the file.
#'  For a given set of arguments to the \code{cluster()} command from within
#'  \emph{mothur}, a number of OTU-clustering results are returned in the same
#'  list file. The exact cutoff values used by \emph{mothur} can vary depending
#'  on the input data. This simple function returns the cutoffs that were actually
#'  included in the \emph{mothur} output. This an important extra step prior to
#'  importing the OTUs with the \code{import_mothur_otulist()} function.
#' 
#' @export
#'  
show_mothur_list_cutoffs <- function(mothur_list_file){
	mothurlist <- readLines(mothur_list_file)
	tabsplit   <- strsplit(mothurlist, "\t", fixed=TRUE)
	cutoffs    <- sapply(tabsplit, function(i){ as.character(i[[1]][1]) })
	return(cutoffs)
}
################################################################################
#' Import mothur list file and return as list object in R.
#'
#' This is a user-available module of a more comprehensive function for importing
#' OTU clustering/abundance data using the \emph{mothur} package. The list object
#' returned by this function is not immediately useable by other \emph{phyloseq}
#' functions, and must be first parsed in conjunction with a separate \emph{mothur}
#' \code{"group"} file. This function is made accessible to \emph{phyloseq} users
#' for troubleshooting and inspection, but the \code{link{import_mothur()}} function
#' is suggested if the goal is to import the OTU clustering results from \emph{mothur}
#' into \emph{phyloseq}.
#'
#' @usage import_mothur_otulist(mothur_list_file, cutoff=NULL)
#'
#' @param mothur_list_file The list file name and/or location as produced by \emph{mothur}.
#'
#' @param cutoff A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur}. The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_list_cutoffs}} function.
#'
#' @return A list, where each element is a character vector of 1 or more 
#'  sequence identifiers, indicating how each sequence from the original data
#'  is clustered into OTUs by \emph{mothur}. Note that in some cases this is highly
#'  dependent on the choice for \code{cutoff}. 
#'
#' @seealso \code{\link{show_mothur_list_cutoffs}}, \code{\link{import_mothur}}
#' 
#' @export
#'  
import_mothur_otulist <- function(mothur_list_file, cutoff=NULL){
	mothurlist <- readLines(mothur_list_file)
	tabsplit   <- strsplit(mothurlist, "\t", fixed=TRUE)
	cutoffs    <- sapply(tabsplit, function(i){ as.character(i[[1]][1]) })
	names(tabsplit) <- cutoffs
	
	# Need to select a cutoff if none was provided by user. 
	# Take the largest non-"unique" cutoff possible,
	# if "unique" is the only cutoff included in the list file, use that.
	if( is.null(cutoff) ){
		if( length(cutoffs) >= 2 ){
			selectCutoffs <- as(cutoffs[cutoffs != "unique"], "numeric")
			cutoff <- as.character(max(selectCutoffs))
		} else {
			# There is only one cutoff value. Use that one.
			cutoff <- cutoffs
		}
	}
	
	# The first two elements are the cutoff and the number of OTUs. metadata. rm.
	OTUs <- tabsplit[[cutoff]][-(1:2)]
	
	# split each element on commas
	OTUs <- strsplit(OTUs, ",", fixed=TRUE)
	
	# Name each OTU, and return the list
	names(OTUs) <- paste("OTUID_", 1:length(OTUs), sep="")
	
	# return as-is
	return(OTUs)
}
################################################################################
#' Parse mothur group file into a simple hash table.
#'
#' The data.frame object
#' returned by this function is not immediately useable by other \emph{phyloseq}
#' functions, and must be first parsed in conjunction with a separate \emph{mothur}
#' \code{"list"} file. This function is made accessible to \emph{phyloseq} users
#' for troubleshooting and inspection, but the \code{link{import_mothur()}} function
#' is suggested if the goal is to import the OTU clustering results from \emph{mothur}
#' into \emph{phyloseq}. You will need both a group file and a list file for that end.
#'
#' @usage import_mothur_groups(mothur_group_file)
#'
#' @param mothur_group_file A character string indicating the location of the 
#'  \emph{mothur}-produced group file in which the sample-source of each sequence
#'  is recorded. See 
#'  \code{http://www.mothur.org/wiki/Make.group}
#' 
#' @return A data.frame that is effectively a hash table between sequence names
#'  and their sample source.
#' 
#' @seealso \code{\link{show_mothur_list_cutoffs}}, \code{\link{import_mothur_otulist}}
#'
#' @export
#'
import_mothur_groups <- function(mothur_group_file){
	group_table <- read.table(mothur_group_file, sep="\t")
	rownames(group_table) <- as(group_table[, 1], "character")
	return(group_table)
}
################################################################################
#' Import mothur list file and return as list object in R.
#'
#' This is a user-available module of a more comprehensive function for importing
#' OTU clustering/abundance data using the \emph{mothur} package. The list object
#' returned by this function is not immediately useable by other \emph{phyloseq}
#' functions, and must be first parsed in conjunction with a separate \emph{mothur}
#' \code{"group"} file. This function is made accessible to \emph{phyloseq} users
#' for troubleshooting and inspection, but the \code{link{import_mothur()}} function
#' is suggested if the goal is to import the OTU clustering results from \emph{mothur}
#' into \emph{phyloseq}.
#'
#' @usage import_mothur_otutable(mothur_list_file, mothur_group_file, cutoff=NULL)
#'
#' @param mothur_list_file The list file name and/or location as produced by \emph{mothur}.
#'
#' @param mothur_group_file The name/location of the group file produced 
#'  by \emph{mothur}'s \code{make.group()} function. It contains information
#'  about the sample source of individual sequences, necessary for creating a
#'  species/taxa abundance table (\code{otuTable}). See
#'  \code{http://www.mothur.org/wiki/Make.group}
#'
#' @param cutoff A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur} (and specified
#'  by the \code{mothur_list_file} argument). The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_list_cutoffs}} function.
#'
#' @return A list, where each element is a character vector of 1 or more 
#'  sequence identifiers, indicating how each sequence from the original data
#'  is clustered into OTUs by \emph{mothur}. Note that in some cases this is highly
#'  dependent on the choice for \code{cutoff}. 
#'
#' @seealso \code{\link{show_mothur_list_cutoffs}}, \code{\link{import_mothur_otulist}}
#' 
#' @export
#'  
import_mothur_otutable <- function(mothur_list_file, mothur_group_file, cutoff=NULL){
	
	otulist       <- import_mothur_otulist(mothur_list_file, cutoff)
	mothur_groups <- import_mothur_groups(mothur_group_file)
	
	# Initialize abundance matrix
	mothur_otu_table <- matrix(0, nrow=length(otulist), ncol=length(levels(mothur_groups[, 2])))
	colnames(mothur_otu_table) <- levels(mothur_groups[, 2])
	rownames(mothur_otu_table) <- names(otulist)
	
	# cycle through each otu, and sum the number of seqs observed for each sample
	for( i in names(otulist) ){
		icount <- tapply(otulist[[i]], mothur_groups[otulist[[i]], 2], length)
		mothur_otu_table[i, names(icount)] <- icount
	}
	
	# Rather than fix the tapply-induced NAs, just replace NAs with 0.
	mothur_otu_table[is.na(mothur_otu_table)] <- 0
	
	# Finally, return the otuTable as a phyloseq otuTable object.
	return(otuTable(mothur_otu_table, speciesAreRows=TRUE))
}
################################################################################
#' Import and prune mothur-produced tree.
#'
#' The \code{ape::read.tree()} function is sufficient for importing a 
#' \emph{mothur}-produced tree into \code{R} in \code{"phylo"} format. This 
#' function further requires a list file as input, and prunes / renames the tree
#' such that it is compatible with the associated abundance data. The expected
#' tree output from \emph{mothur} will contain a tip for each sequence in the
#' original dataset, even if those sequences are highly similar or identical. 
#' This function returns a reduced tree that has been pruned to contain only unique
#' taxa, as defined by a particular OTU-clustering cutoff in the \emph{mothur}
#' list file.
#'
#' This is a user-available module of a more comprehensive function for importing
#' output files from the \emph{mothur} package, \code{link{import_mothur}}. 
#' The \code{link{import_mothur}} function
#' is suggested if the goal is to import more than just the tree from \emph{mothur}.
#'
#' @usage import_mothur_tree(mothur_tree_file, mothur_list_file, cutoff=NULL)
#'
#' @param mothur_tree_file The tree file name that was output from \emph{mothur}.
#'  Probably a file that ends with the suffix \code{".tree"}.
#'
#' @param mothur_list_file The list file name and/or location as produced by \emph{mothur}.
#'
#' @param cutoff A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur} (and specified
#'  by the \code{mothur_list_file} argument). The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_list_cutoffs}} function.
#'
#' @return A reduced tree that has been pruned to contain only unique
#'  taxa, as defined by a particular OTU-clustering cutoff in the \emph{mothur}
#'  list file. 
#'
#' @seealso \code{\link{import_mothur}}, \code{\link{import_mothur_otulist}}
#' 
#' @export
#'
import_mothur_tree <- function(mothur_tree_file, mothur_list_file, cutoff=NULL){
	
	otulist <- import_mothur_otulist(mothur_list_file, cutoff)

	# Read the original all-sequences tree from mothur
	# # original_tree <- read.tree(mothur_tree_file)
	# # tree          <- original_tree
	tree <- read.tree(mothur_tree_file)
	
	# Loop to merge the sequences in the tree by OTU
	# cycle through each otu, and sum the number of seqs observed for each sample
	for( i in names(otulist) ){
		# i <- names(otulist)[1]
		# First merge the reads that are in the same OTU ("eqspecies" argument)
		tree <- mergespecies(tree, otulist[[i]])
		# Rename the tip that was kept to the otuID.
		# By default, the first element of eqspecies is used as archetype.
		# This also ensures reliable behavior in the instances of singleton OTUs
		tree$tip.label[tree$tip.label == otulist[[i]][1]] <- i
	}
	return(tree)
}
################################################################################
#' General function for importing mothur files into phyloseq.
#'
#' @usage import_mothur(mothur_list_file,  mothur_group_file=NULL, mothur_tree_file=NULL, cutoff=NULL)
#'
#' @param mothur_list_file Required. The list file name / location produced by \emph{mothur}.
#'
#' @param mothur_group_file Optional. The name/location of the group file produced 
#'  by \emph{mothur}'s \code{make.group()} function. It contains information
#'  about the sample source of individual sequences, necessary for creating a
#'  species/taxa abundance table (\code{otuTable}). See
#'  \code{http://www.mothur.org/wiki/Make.group} 
#'
#' @param mothur_tree_file Optional. The tree file name produced by \emph{mothur}.
#'  Probably a file that ends with the suffix \code{".tree"}.
#'
#' @param cutoff A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur} (and specified
#'  by the \code{mothur_list_file} argument). The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_list_cutoffs}} function.
#'
#' @return The object class depends on the provided arguments.
#'  If the first three arguments are provided, then an \code{otuTree} object should
#'  be returned, containing both an OTU-only tree and its associated 
#'  \code{otuTable}-class abundance table. If only a list and group file are 
#'  provided, then an \code{otuTable} object is returned. Similarly, if only a list
#'  and tree file are provided, then only a tree is returned (\code{"phylo"} class).
#'
#' @export
#'
#' @examples
#' # # The following example assumes you have downloaded the esophagus example
#' # # dataset from the mothur wiki:
#' # # "http://www.mothur.org/wiki/Esophageal_community_analysis"
#' # # "http://www.mothur.org/w/images/5/55/Esophagus.zip"
#' # # The path on your machine may (probably will) vary
#' # mothur_list_file  <- "~/Downloads/mothur/Esophagus/esophagus.an.list"
#' # mothur_group_file <- "~/Downloads/mothur/Esophagus/esophagus.good.groups"
#' # mothur_tree_file  <- "~/Downloads/mothur/Esophagus/esophagus.tree"
#' # # # Actual examples follow:
#' # show_mothur_list_cutoffs(mothur_list_file)
#' # test1 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file)
#' # test2 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file, cutoff="0.02")
#' # # Returns just a tree
#' # import_mothur(mothur_list_file, mothur_tree_file=mothur_tree_file)
#' # # Returns just an otuTable
#' # import_mothur(mothur_list_file, mothur_group_file=mothur_group_file)
#' # # Returns an error
#' # import_mothur(mothur_list_file)
#' # # Should return an "OMG, you must provide the list file" error
#' # import_mothur()
import_mothur <- function(mothur_list_file,  mothur_group_file=NULL,
		mothur_tree_file=NULL, cutoff=NULL){

	if( missing(mothur_list_file) ){cat("you must provide the mothur_list_file argument\n")}

	# Create otuTable object, OTU, only if group file provided. 
	if( !is.null(mothur_group_file) ){
		OTU  <- import_mothur_otutable(mothur_list_file, mothur_group_file, cutoff)		
	} 
	
	# Similarly, only get modified tree object if tree file provided.
	if( !is.null(mothur_tree_file) ){
		tree <- import_mothur_tree(mothur_tree_file, mothur_list_file, cutoff)		
	}

	### Decide what is going to be returned, based on arguments.
	# If no group and no tree, then only list. Return list module output.
	if( is.null(mothur_group_file) & is.null(mothur_tree_file) ){
		return(import_mothur_otulist(mothur_list_file, cutoff))
	} else if( !is.null(mothur_group_file) & is.null(mothur_tree_file) ){
		return(OTU)
	} else if( is.null(mothur_group_file) & !is.null(mothur_tree_file) ){
		return(tree)
	} else if( !is.null(mothur_group_file) & !is.null(mothur_tree_file) ){
		# Return a merged OTU and tree object
		return(phyloseq(OTU, tree))		
	}
}
################################################################################
################################################################################