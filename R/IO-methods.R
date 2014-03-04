################################################################################
#' Universal import method (wrapper) for phyloseq-package
#'
#' A user must still understand the additional arguments required for each
#' type of import data. Those arguments are described in detail at the
#' tool-specific \code{import_*} links below. Each clustering tool / package / pipeline
#' has its own idiosyncratic set of file names / types, and it remains the
#' responsibility of the user to understand which file-path should be provided
#' to each argument for the particular importing submethod. This method
#' merely provides a central documentation and method-name, and the arguments
#' are passed along as-is.
#'
#' @usage import(pipelineName, ...)
#'
#' @param pipelineName (Required). Character string. The name of the 
#'  analysis tool / pipeline / package
#'  that created the OTU-cluster data or other data that you now want to import.
#'  Current options are \code{c("mothur", "pyrotagger", "QIIME", "RDP")}, and
#'  only the first letter is necessary.
#'
#' @param ... (Required). Additional named arguments providing file paths, and possible
#'  other paramaters to the desired tool-specific import function.
#'
#' @return In most cases a \code{\link{phyloseq-class}} will be returned, though
#'  the included component data will vary by pipeline/tool, and also
#'  by the types of data files provided.
#'  The expected behavior is to return the most-comprehensive object possible,
#'  given the provided arguments and pipeline/tool.
#'
#' @seealso 
#' 
#' For BIOM format, see:
#' \code{\link{import_biom}}
#' 
#' For mothur, see:
#' \code{\link{import_mothur}}
#' 
#' Separate tools for mothur are also:
#' \code{\link{show_mothur_cutoffs}}
#' \code{\link{import_mothur_dist}}
#' \code{\link{export_mothur_dist}}
#' 
#' For PyroTagger, see:
#' \code{\link{import_pyrotagger_tab}}
#'
#' For QIIME legacy format, see:
#' \code{\link{import_qiime}}
#' 
#' For RDP pipeline, see:
#' \code{\link{import_RDP_cluster}}
#'
#' \code{\link{import_RDP_otu}}
#' 
#' @references  
#'
#' BIOM: \url{http://www.biom-format.org/}
#' 
#' mothur: \url{http://www.mothur.org/wiki/Main_Page}
#'
#' PyroTagger: \url{http://pyrotagger.jgi-psf.org/}
#'
#' QIIME: \url{http://qiime.org/}
#'
#' RDP pipeline: \url{http://pyro.cme.msu.edu/index.jsp}
#' 
#' @export
#' @examples
#'  ## See documentation of a specific import function
import <- function(pipelineName, ...){
	# Reduce pipelineName to just its first letter, as all are different
	pipelineName <- substr(pipelineName, 1, 1)

	# Test that it is in the set
	if( !(pipelineName %in% c("B", "b", "M", "m", "P", "p", "Q", "q", "R", "r")) ){
		stop("You need to select among available importer types:\n", 
		"\"BIOM\", \"mothur\", \"pyrotagger\", \"QIIME\", \"RDP\" \n See ?import for details")
	}

	if( pipelineName %in% c("B", "b") ){
		return( import_biom(...) )
	}	
	if( pipelineName %in% c("M", "m") ){
		return( import_mothur(...) )
	}
	if( pipelineName %in% c("P", "p") ){
		return( import_pyrotagger_tab(...) ) 
	}	
	if( pipelineName %in% c("Q", "q") ){
		return( import_qiime(...) )
	}
	if( pipelineName %in% c("R", "r") ){
		return( import_RDP_cluster(...) )
	}
}
################################################################################
#' Import function to read the now legacy-format QIIME OTU table.
#'
#' QIIME produces several files that can be directly imported by
#' the \code{\link{phyloseq-package}}.
#' Originally, QIIME produced its own custom format table
#' that contained both OTU-abundance
#' and taxonomic identity information.
#' This function is still included in phyloseq mainly to accommodate these
#' now-outdated files. Recent versions of QIIME store output in the
#' biom-format, an emerging file format standard for microbiome data.
#' If your data is in the biom-format, if it ends with a \code{.biom}
#' file name extension, then you should use the \code{\link{import_biom}}
#' function instead.
#' 
#' Other related files include
#' the mapping-file that typically stores sample covariates,
#' converted naturally to the 
#' \code{\link{sample_data-class}} component data type in the phyloseq-package. 
#' QIIME may also produce a
#' phylogenetic tree with a tip for each OTU, which can also be imported
#' specified here or imported separately using \code{\link{read_tree}}.
#' 
#' See \url{"http://www.qiime.org/"} for details on using QIIME. While there are
#' many complex dependencies, QIIME can be downloaded as a pre-installed 
#' linux virtual machine that runs ``off the shelf''. 
#'
#' The different files useful for import to \emph{phyloseq} are not collocated in
#' a typical run of the QIIME pipeline. See the main \emph{phyloseq} vignette for an 
#' example of where ot find the relevant files in the output directory. 
#'
#' @param otufilename (Optional). A character string indicating 
#'  the file location of the OTU file.
#'  The combined OTU abundance and taxonomic identification file, 
#'  tab-delimited, as produced by QIIME under default output settings.
#'  Default value is \code{NULL}. 
#' 
#' @param mapfilename (Optional). The QIIME map file is required
#'  for processing barcoded primers in QIIME
#'  as well as some of the post-clustering analysis. This is a required
#'  input file for running QIIME. Its strict formatting specification should be
#'  followed for correct parsing by this function.
#'  Default value is \code{NULL}. 
#'
#' @param treefilename (Optional). Default value is \code{NULL}. 
#'  A file representing a phylogenetic tree
#'  or a \code{\link{phylo}} object.
#'  Files can be NEXUS or Newick format.
#'  See \code{\link{read_tree}} for more details.
#'  Also, if using a recent release of the GreenGenes database tree,
#'  try the \code{\link{read_tree_greengenes}} function --
#'  this should solve some issues specific to importing that tree. 
#'  If provided, the tree should have the same OTUs/tip-labels
#'  as the OTUs in the other files. 
#'  Any taxa or samples missing in one of the files is removed from all. 
#'  As an example from the QIIME pipeline,
#'  this tree would be a tree of the representative 16S rRNA sequences from each OTU
#'  cluster, with the number of leaves/tips equal to the number of taxa/species/OTUs,
#'  or the complete reference database tree that contains the OTU identifiers
#'  of every OTU in your abundance table.
#'  Note that this argument can be a tree object (\code{\link[ape]{phylo}}-class)
#'  for cases where the tree has been --- or needs to be --- imported separately,
#'  as in the case of the GreenGenes tree mentioned earlier (code{\link{read_tree_greengenes}}).
#'
#' @param refseqfilename (Optional). Default \code{NULL}.
#'  The file path of the biological sequence file that contains at a minimum
#'  a sequence for each OTU in the dataset.
#'  Alternatively, you may provide an already-imported
#'  \code{\link[Biostrings]{XStringSet}} object that satisfies this condition.
#'  In either case, the \code{\link{names}} of each OTU need to match exactly the
#'  \code{\link{taxa_names}} of the other components of your data.
#'  If this is not the case, for example if the data file is a FASTA format but
#'  contains additional information after the OTU name in each sequence header,
#'  then some additional parsing is necessary,
#'  which you can either perform separately before calling this function, 
#'  or describe explicitly in a custom function provided in the (next) argument,
#'  \code{refseqFunction}.
#'  Note that the \code{\link[Biostrings]{XStringSet}} class can represent any 
#'  arbitrary sequence, including user-defined subclasses, but is most-often
#'  used to represent RNA, DNA, or amino acid sequences. 
#'  The only constraint is that this special list of sequences
#'  has exactly one named element for each OTU in the dataset.
#' 
#' @param refseqFunction (Optional). 
#'  Default is \code{\link[Biostrings]{readDNAStringSet}}, 
#'  which expects to read a fasta-formatted DNA sequence file.
#'  If your reference sequences for each OTU are amino acid, RNA, or something else,
#'  then you will need to specify a different function here.
#'  This is the function used to read the file connection provided as the
#'  the previous argument, \code{refseqfilename}.
#'  This argument is ignored if \code{refseqfilename} is already a
#'  \code{\link[Biostrings]{XStringSet}} class. 
#' 
#' @param refseqArgs (Optional). 
#'  Default \code{NULL}.
#'  Additional arguments to \code{refseqFunction}.
#'  See \code{\link[Biostrings]{XStringSet-io}} for details about
#'  additional arguments to the standard read functions in the Biostrings package.
#'
#' @param parseFunction (Optional). An optional custom function for parsing the
#'  character string that contains the taxonomic assignment of each OTU. 
#'  The default parsing function is \code{\link{parse_taxonomy_qiime}},
#'  specialized for splitting the \code{";"}-delimited strings and also 
#'  attempting to interpret greengenes prefixes, if any, as that is a common
#'  format of the taxonomy string produced by QIIME.
#'  
#' @param verbose (Optional). A \code{\link{logical}}.
#'  Default is \code{TRUE}. 
#'  Should progresss messages
#'  be \code{\link{cat}}ted to standard out?
#'
#' @param ... Additional arguments passed to \code{\link{read_tree}}
#'
#' @return A \code{\link{phyloseq-class}} object.
#'
#' @seealso
#'
#' \code{\link{phyloseq}}
#'
#' \code{\link{merge_phyloseq}}
#' 
#' \code{\link{read_tree}}
#' 
#' \code{\link{read_tree_greengenes}}
#'
#' \code{\link[Biostrings]{XStringSet-io}}
#'
#' @references \url{http://qiime.org/}
#'
#' ``QIIME allows analysis of high-throughput community sequencing data.''
#' J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman,
#' Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon,
#' Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, 
#' Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky,
#' Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight;
#' Nature Methods, 2010; doi:10.1038/nmeth.f.303
#'
#' @importClassesFrom Biostrings XStringSet
#' @importFrom Biostrings readDNAStringSet
#' @export
#' @examples
#'  otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
#'  mapfile <- system.file("extdata", "master_map.txt", package="phyloseq")
#'  trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
#'  import_qiime(otufile, mapfile, trefile)
import_qiime <- function(otufilename=NULL, mapfilename=NULL,
	treefilename=NULL, refseqfilename=NULL, 
  refseqFunction=readDNAStringSet, refseqArgs=NULL,
	parseFunction=parse_taxonomy_qiime, verbose=TRUE, ...){

	# initialize the argument-list for phyloseq. Start empty.
	argumentlist <- list()

	if( !is.null(mapfilename) ){	
		if( verbose ){
			cat("Processing map file...", fill=TRUE)
		}
		QiimeMap     <- import_qiime_sample_data(mapfilename)
		argumentlist <- c(argumentlist, list(QiimeMap))
	}

	if( !is.null(otufilename) ){
		if( verbose ){
			cat("Processing otu/tax file...", fill=TRUE)
		}		
		otutax <- import_qiime_otu_tax(otufilename, parseFunction, verbose=verbose)
		otutab <- otu_table(otutax$otutab, TRUE)
		taxtab <- tax_table(otutax$taxtab)
		argumentlist <- c(argumentlist, list(otutab), list(taxtab) )
	}

	if( !is.null(treefilename) ){
		if(verbose){cat("Processing phylogenetic tree...\n", treefilename, "...\n")}
		if(inherits(treefilename, "phylo")){
			# If argument is already a tree, don't read, just assign.
			tree = treefilename
		} else {
      # If it is not a tree, assume file and attempt to import.
			# NULL is silently returned if tree is not read properly.
			tree <- read_tree(treefilename, ...)
		}
		# Add to argument list or warn
		if( is.null(tree) ){
			warning("treefilename failed import. It will not be included.")
		} else {
			argumentlist <- c(argumentlist, list(tree) )
		}
	}
	
	if( !is.null(refseqfilename) ){
		if( verbose ){
			cat("Processing Reference Sequences...", fill=TRUE)
		}
		if( inherits(refseqfilename, "XStringSet") ){
			# If argument is already a XStringSet, don't read, just assign.
			refseq = refseqfilename
		} else {
			# call refseqFunction and read refseqfilename, 
      # either with or without additional args
			if( !is.null(refseqArgs) ){
				refseq = do.call("refseqFunction", c(list(refseqfilename), refseqArgs))
			} else {
				refseq = refseqFunction(refseqfilename)
			}
		}
		argumentlist <- c(argumentlist, list(refseq) )
	}

	do.call("phyloseq", argumentlist)
}
################################################################################
#' Somewhat flexible tree-import function
#'
#' This function is a convenience wrapper around the
#' \code{\link[ape]{read.tree}} (Newick-format) and
#' \code{\link[ape]{read.nexus}} (Nexus-format) importers provided by
#' the \code{\link[ape]{ape-package}}. This function attempts to return a valid
#' tree if possible using either format importer. If it fails, it silently 
#' returns \code{NULL} by default, rather than throwing a show-stopping error.
#'
#' @usage read_tree(treefile, errorIfNULL=FALSE, ...)
#'
#' @param treefile (Required). A character string implying a file \code{\link{connection}}
#'  (like a path or URL), or an actual \code{\link{connection}}.
#'  Must be a Newick- or Nexus-formatted tree.
#'
#' @param errorIfNULL (Optional). Logical. Should an error be thrown if no tree
#'  can be extracted from the connection?
#'  Default is \code{FALSE}, indicating that \code{NULL} will be 
#'  SILENTLY returned, rather than an error. 
#'  Be cautious about this behavior. Useful for phyloseq internals, but might
#'  be hard to track in your own code if you're not aware of this
#'  ``no error by default'' setting. If this is a problem, change this value
#'  to \code{TRUE}, and you can still use the function.
#'
#' @param ... (Optional). Additional parameter(s) passed to the
#'  relevant tree-importing function.
#'
#' @return If successful, returns a \code{\link{phylo}}-class object as defined
#'  in the \code{\link[ape]{ape-package}}. Returns NULL if neither tree-reading function worked.
#'
#' @seealso
#' \code{\link{read_tree_greengenes}}
#'
#' \code{\link{phylo}}
#' 
#' \code{\link[ape]{read.tree}}
#'  
#' \code{\link[ape]{read.nexus}}
#'
#' @importFrom ape read.nexus
#' @importFrom ape read.tree
#' @export
#' @examples
#' read_tree(system.file("extdata", "esophagus.tree.gz", package="phyloseq"))
#' read_tree(system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq"))
read_tree <- function(treefile, errorIfNULL=FALSE, ...){
	# "phylo" object provided directly
	if( class(treefile)[1] %in% c("phylo") ){ 
		tree <- treefile
	} else {
	  # file path to tree file provided.
		# Try Nexus first, protected, then newick if it fails
		tree <- NULL
		try(tree <- read.nexus(treefile, ...), TRUE)
		# Try Newick if nexus didn't work.
		if(is.null(tree)) try(tree <- read.tree(treefile, ...), TRUE)
	}
	# If neither tree-import worked (still NULL), report warning
	if( errorIfNULL & is.null(tree) ){
		stop("tree file could not be read.\nPlease retry with valid tree.")
	}
  if( !is.null(tree) ){
    # Perform any standard phyloseq checks/fixes
    # E.g. Replace any NA branch-length values in the tree with zero.
    tree = fix_phylo(tree)    
  }
	return(tree)
}
################################################################################
#' Read GreenGenes tree released in annotated newick format
#'
#' In principal, this is a standard newick format, that can be imported
#' into R using \code{\link{read_tree}},
#' which in-turn utilizes \code{\link[ape]{read.tree}}.
#' However, \code{\link[ape]{read.tree}} has failed to import
#' recent (October 2012 and later) releases of the GreenGenes tree,
#' and this problem has been traced to the additional annotations
#' added to some internal nodes 
#' that specify taxonomic classification between single-quotes. 
#' To solve this problem and create a clear container
#' for fixing future problems with the format of GreenGenes-released trees,
#' this function is available in phyloseq and exported for users.
#' It is also referenced in the documentation of the import functions
#' for QIIME legacy and BIOM format importers --
#' \code{\link{import_qiime}} and \code{\link{import_biom}}, respectively.
#' However, since the precise format of the tree is not restricted to GreenGenes trees
#' by QIIME or for the biom-format, this function is not called
#' automatically by those aforementioned import functions.
#' If your tree is formatted like, or is one of, the official GreenGenes
#' release trees, then you should use this function and provide its output
#' to your relevant import function.
#' 
#' @param treefile (Required). A character string implying 
#'  a file \code{\link{connection}}
#'  (like a path or URL), or an actual \code{\link{connection}}.
#'  Must be a Newick--formatted tree released by GreenGenes
#'  in October 2012 or later.
#'  The similarity threshold of the OTUs should not matter,
#'  except that it should match your OTU table.
#'  
#' @return A tree, represented as a \code{\link{phylo}} object.
#' 
#' @importFrom ape read.tree
#' @export
#' @examples
#' # Read the May 2013, 73% similarity official tree,
#' # included as extra data in phyloseq.
#' treefile = system.file("extdata", "gg13-5-73.tree.gz", package="phyloseq")
#' x = read_tree_greengenes(treefile)
#' x
#' class(x)
#' y = read_tree(treefile)
#' y
#' class(y)
#' ## Not run, causes an error:
#' # library("ape")
#' # read.tree(treefile)
read_tree_greengenes = function(treefile){
  alines = readLines(treefile, warn=FALSE)
  # Collapse to one line, in case it isn't already.
  alines = paste0(alines, collapse="")
  # replace all semicolons with something weird
  # that isn't already a special newick character.
  newdelim = "><-><"
  clines = gsub("\\;", newdelim, alines)
  # reinstate the final character as a semicolon
  clines = gsub(paste0(newdelim, "$"), ";", clines)
  # Convert your newick string into a phylo-class tree.
  tree = read.tree("", text=clines)
  # Now that it is phylo-class, reinstate semicolon
  # as the delimiter in the node labels
  gsub(newdelim, ";", tree$node.label)
  # Also get rid of those extra quotes
  gsub("'", "", tree$node.label)
  # Return the cleaned-up tree
  return(tree)
}
################################################################################
#' Import now legacy-format QIIME OTU table as a list of two matrices.
#'
#' Now a legacy-format, older versions of QIIME
#' produced an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information in a tab-delimted table.
#' If your file ends with the extension \code{.biom}, or if you happen to know
#' that it is a biom-format file, or if you used default settings in a version
#' of QIIME of \code{1.7} or greater,
#' then YOU SHOULD USE THE BIOM-IMPORT FUNCTION instead, 
#' \code{\link{import_biom}}.
#'
#' This function uses chunking to perform both the reading and parsing in blocks 
#' of optional size,
#' thus constrain the peak memory usage. 
#' feature should make this
#' importer accessible to machines with modest memory,
#' but with the caveat that
#' the full numeric matrix must be a manageable size at the end, too.
#' In principle, the final tables will be large, but much more efficiently represented than
#' the character-stored numbers.
#' If total memory for storing the numeric matrix becomes problematic,
#' a switch to a sparse matrix representation of the abundance
#' -- which is typically well-suited to this data -- might provide a solution.
#'
#' @param file (Required). The path to the qiime-formatted file you want to
#'  import into R. Can be compressed (e.g. \code{.gz}, etc.), though the
#'  details may be OS-specific. That is, Windows-beware.
#'
#' @param parseFunction (Optional). An optional custom function for parsing the
#'  character string that contains the taxonomic assignment of each OTU. 
#'  The default parsing function is \code{\link{parse_taxonomy_qiime}},
#'  specialized for splitting the \code{";"}-delimited strings and also 
#'  attempting to interpret greengenes prefixes, if any, as that is a common
#'  format of the taxonomy string produced by QIIME.
#'
#' @param verbose (Optional). A \code{\link{logical}}.
#'  Default is \code{TRUE}. 
#'  Should progresss messages
#'  be \code{\link{cat}}ted to standard out?
#'
#' @param parallel (Optional). Logical. Should the parsing be performed in 
#'  parallel?. Default is \code{FALSE}. Only a few steps are actually 
#'  parallelized, and for most datasets it will actually be faster and 
#'  more efficient to keep this set to \code{FALSE}.
#'  Also, to get any benefit at all, you will need to register a 
#'  parallel ``backend'' through one of the backend packages supported
#'  by the \code{\link{foreach-package}}.
#'  
#' @return A list of two matrices. \code{$otutab} contains the OTU Table
#'  as a numeric matrix, while \code{$taxtab} contains a character matrix
#'  of the taxonomy assignments.
#'
#' @importFrom data.table fread
#' @importFrom plyr llply
#'
#' @seealso
#' \code{\link{import}}
#' 
#' \code{\link{merge_phyloseq}}
#' 
#' \code{\link{phyloseq}}
#' 
#' \code{\link{import_qiime}}
#' 
#' \code{\link{read_tree}}
#' 
#' \code{\link{read_tree_greengenes}}
#' 
#' \code{\link{import_env_file}}
#' 
#' @export
#' @examples
#'  otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
#'  import_qiime_otu_tax(otufile)
import_qiime_otu_tax <- function(file, parseFunction=parse_taxonomy_qiime,
                                 verbose=TRUE, parallel=FALSE){
  if(verbose){cat("Reading file into memory prior to parsing...\n")}
  x = readLines(file)
  if(verbose){cat("Detecting first header line...\n")}
  # Check for commented lines, starting with line 1.
  # The deepest commented line (largest n) is assumed to have header information.
  skipLines = max(which(substr(x[1:25L], 1, 1)=="#"))-1L
  if(verbose){cat("Header is on line", (skipLines + 1L), " \n")}
  if(verbose){cat("Converting input file to a table...\n")}
  x = fread(input=paste0(x, collapse="\n"), sep="\t", header=TRUE, skip=skipLines)
  if(verbose){cat("Defining OTU table... \n")}
  taxstring = x$`Consensus Lineage`
  # This pops the taxonomy (Consensus Lineage) column, in-place statement
  x[, `Consensus Lineage`:=NULL]
  # Store the OTU names, you will pop the column
  OTUnames = x$`#OTU ID`
  # This pops the OTUID column, in-place statement
  x[, `#OTU ID`:=NULL]
  x <- as(x, "matrix")
  rownames(x) <- OTUnames
  rm(OTUnames)
  if(verbose){cat("Parsing taxonomy table...\n")}
  # Split into "jagged" list (vectors of different lengths)
  taxlist = llply(taxstring, parseFunction, .parallel=parallel)
  # Add OTU names to list element names
  names(taxlist) <- rownames(x)
  # Build the tax table from the jagged list.  
  taxtab <- build_tax_table(taxlist)
	# Call garbage collection one more time. Lots of unneeded stuff.
	garbage.collection <- gc(FALSE)
  # Return the named list
	return(list(otutab=x, taxtab=taxtab))
}
################################################################################
################################################################################
#' Import just \code{sample_data} file from QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' This includes the map-file, which is an important \emph{input}
#' to QIIME that can also indicate sample covariates. It is converted naturally to the 
#' sample_data component data type in phyloseq-package, based on the R data.frame.
#' 
#' See \code{\link{import_qiime}} for more information about QIIME. It is also the
#' suggested function for importing QIIME-produced data files. 
#'
#' @usage import_qiime_sample_data(mapfilename)
#'
#' @param mapfilename (Required). A character string or connection.
#'  That is, any suitable \code{file} argument to the \code{\link{read.table}} function. 
#'  The name of the QIIME map
#'  file required for processing pyrosequencing tags
#'  in QIIME as well as some of the post-clustering analysis. This is a required
#'  input file for running QIIME. Its strict formatting specification is expected by
#'  this function, do not attempt to modify it manually once it has worked properly
#'  in QIIME. 
#'
#' @return A \code{sample_data} object.
#'
#' @seealso
#' 
#' \code{\link{import}}
#' 
#' \code{\link{merge_phyloseq}}
#' 
#' \code{\link{phyloseq}}
#' 
#' \code{\link{import_qiime}}
#' 
#' \code{\link{import_qiime_otu_tax}}
#' 
#' \code{\link{import_env_file}}
#' 
#' @export
#' @examples 
#'  mapfile <- system.file("extdata", "master_map.txt", package = "phyloseq")
#'  import_qiime_sample_data(mapfile)
import_qiime_sample_data <- function(mapfilename){
	# Process mapfile. Name rows as samples.
	QiimeMap <- read.table(file=mapfilename, header=TRUE,
		sep="\t", comment.char="")
	rownames(QiimeMap) <- as.character(QiimeMap[,1])
	return( sample_data(QiimeMap) )
}
################################################################################
#' Read a UniFrac-formatted ENV file.
#'
#' Convenience wrapper function to read the environment-file, as formatted for 
#' input to the UniFrac server (\url{http://bmf2.colorado.edu/unifrac/}).
#' The official format of these files is that
#' each row specifies (in order) the sequence name, source sample, and (optionally)
#' the number of times the sequence was observed. 
#'
#' @usage import_env_file(envfilename, tree=NULL, sep="\t", ...)
#'
#' @param envfilename (Required). A charater string of the ENV filename (relative or absolute)
#'
#' @param tree (Optional). \code{\link{phylo-class}} object to be paired with
#'  the output otu_table. 
#'
#' @param sep A character string indicating the delimiter used in the file.
#'  The default is \code{"\t"}.
#'
#' @param ... Additional parameters passed on to \code{\link{read.table}}.
#'
#' @return An \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}} if 
#'  a \code{\link{phylo-class}} argument is provided to \code{tree}.
#'
#' @references \url{http://bmf2.colorado.edu/unifrac/}
#' 
#' @seealso
#' \code{\link{import}}
#' 
#' \code{\link{tip_glom}}
#' @export
#' @examples 
#' # import_env_file(myEnvFile, myTree)
import_env_file <- function(envfilename, tree=NULL, sep="\t", ...){
	tipSampleTable <- read.table(envfilename, sep=sep, ...)
	# Convert to otu_table-class table (trivial table)
	physeq <- envHash2otu_table(tipSampleTable)
	# If tree is provided, combine it with the OTU Table
	if( class(tree) == "phylo" ){
		# Create phyloseq-class with a tree and OTU Table (will perform any needed pruning)
		physeq <- phyloseq(physeq, tree)
	}
	return(physeq)
}
################################################################################  
#' Convert a sequence-sample hash (like ENV file) into an OTU table.
#' 
#' Parses an ENV-file into a sparse matrix of species-by-sample, where
#' each species-row has only one non-zero value. We call this sparse abundance
#' table the trivial OTU table, where every sequence is treated as a separate 
#' species. If a phylogenetic tree is available, it can be submitted with this
#' table as arguments to \code{\link{tip_glom}} to create an object with a
#' non-trivial \code{otu_table}.  
#'
#' @usage envHash2otu_table(tipSampleTable)
#'
#' @param tipSampleTable (Required). A two-column character table (matrix or data.frame), 
#' where each row specifies the sequence name and source sample, consistent with the 
#' env-file for the UniFrac server (\url{http://bmf2.colorado.edu/unifrac/}). 
#'
#' @return \code{\link{otu_table}}. A trivial OTU table where each sequence 
#'  is treated as a separate OTU. 
#' 
#' @references \url{http://bmf2.colorado.edu/unifrac/}
#' 
#' @seealso
#' \code{\link{import_env_file}}
#' 
#' \code{\link{tip_glom}}
#' 
#' \code{\link{otu_table}}
#'
#' @keywords internal
#' @examples #
#' ## fakeSeqNameVec <- paste("seq_", 1:8, sep="")
#' ## fakeSamNameVec <- c(rep("A", 4), rep("B", 4))
#' ## fakeSeqAbunVec <- sample(1:50, 8, TRUE)
#' ## test    <- cbind(fakeSeqNameVec, fakeSamNameVec, fakeSeqAbunVec)
#' ## testotu <- envHash2otu_table( test )
#' ## test    <- cbind(fakeSeqNameVec, fakeSamNameVec)
#' ## testotu <- envHash2otu_table( test )
envHash2otu_table <- function(tipSampleTable){
	if( ncol(tipSampleTable) > 2 ){
		tst <- tipSampleTable
		trivialOTU <- matrix(0, nrow=nrow(tst), ncol=length(unique(tst[,2])))
		colnames(trivialOTU) <- unique(tst[,2])
		rownames(trivialOTU) <- tst[,1]
		for( i in 1:nrow(tst) ){
			trivialOTU[tst[i, 1], tst[i, 2]] <- as.integer(tst[i, 3])
		}
	} else {
		trivialOTU <- table(as.data.frame(tipSampleTable))
		trivialOTU <- as(trivialOTU, "matrix")	
	}
	return( otu_table(trivialOTU, taxa_are_rows=TRUE) )
}
################################################################################
################################################################################
#' Import RDP cluster file and return otu_table (abundance table).
#'
#' The RDP cluster pipeline (specifically, the output of the complete linkage clustering step)
#' has no formal documentation for the \code{".clust"} 
#' file or its apparent sequence naming convention. 
#'
#' \code{http://pyro.cme.msu.edu/index.jsp}
#'
#' The cluster file itself contains
#' the names of all sequences contained in input alignment. If the upstream
#' barcode and aligment processing steps are also done with the RDP pipeline,
#' then the sequence names follow a predictable naming convention wherein each
#' sequence is named by its sample and sequence ID, separated by a \code{"_"} as
#' delimiter:
#'
#' \code{"sampleName_sequenceIDnumber"}
#'
#' This import function assumes that the sequence names in the cluster file follow
#' this convention, and that the sample name does not contain any \code{"_"}. It
#' is unlikely to work if this is not the case. It is likely to work if you used
#' the upstream steps in the RDP pipeline to process your raw (barcoded, untrimmed)
#' fasta/fastq data. 
#'
#' This function first loops through the \code{".clust"} file and collects all
#' of the sample names that appear. It secondly loops through each OTU (\code{"cluster"};
#' each row of the cluster file) and sums the number of sequences (reads) from
#' each sample. The resulting abundance table of OTU-by-sample is trivially
#' coerced to an \code{\link{otu_table}} object, and returned.
#'
#' @usage import_RDP_cluster(RDP_cluster_file)
#'
#' @param RDP_cluster_file A character string. The name of the \code{".clust"} 
#'  file produced by the 
#'  the complete linkage clustering step of the RDP pipeline.
#'
#' @return An \code{\link{otu_table}} object parsed from the \code{".clust"} file. 
#' 
#' @references \url{http://pyro.cme.msu.edu/index.jsp}
#'
#' @export
#' 
import_RDP_cluster <- function(RDP_cluster_file){
	
	# Read file and pop the header lines
	RDP_raw_otu_lines_only <- readLines(RDP_cluster_file)[-(1:5)]
	
	# internal function:
	make_verbose_sample_list <- function(RDP_raw_otu_lines_only){
		# Each OTU line has a 3 element "line header" that indicates the OTUID, the name of the file,
		# and the number of sequences that are included in this cluster.
		# From each line, remove the header elements
		get_sample_names_from_one_line <- function(otuline){
			# first split the line on tabs "\t"
			splittabs <- strsplit(otuline, "\t")[[1]]
			
			# next, remove the header by keeping on the 4th element.
			seqIDonly <- splittabs[4]
			
			# Finally, split on white space
			seqIDonly <- strsplit(seqIDonly, "[[:space:]]+")[[1]]
			
			# For each element in seqIDonly, split on the underscore delimiter
			splitseqnames <- strsplit(seqIDonly, "_", fixed=TRUE)
			
			# Return the sample names from the first element (assumes no "_" in sample names)
			return( sapply(splitseqnames, function(i){i[1]}) )
		}
		return( sapply(RDP_raw_otu_lines_only, get_sample_names_from_one_line) )
	}
	
	## Get the verbose sample name list, and then shrink to the 
	## unique sample names in the entire dataset. 
	## Need this unique list for initializing the OTU abundance matrix
	RDPsamplenameslist <- make_verbose_sample_list(RDP_raw_otu_lines_only)
	RDPsamplenames     <- unique(unlist(RDPsamplenameslist))
	
	# remove NAs
	RDPsamplenames <- RDPsamplenames[!is.na(RDPsamplenames)]
	
	# Initialize otu abundance matrix.
	otumat <- matrix(0, nrow=length(RDP_raw_otu_lines_only), ncol=length(RDPsamplenames))
	rownames(otumat) <- paste("OTUID_", 1:length(RDP_raw_otu_lines_only))
	colnames(otumat) <- RDPsamplenames
	
	# Now re-loop through the cluster file (by OTU) and sum the 
	# abundance of sequences from each sample 
	for( i in 1:length(RDP_raw_otu_lines_only) ){
		# i = 1
	
		# first split the line on tabs "\t"
		splittabs <- strsplit(RDP_raw_otu_lines_only[i], "\t")[[1]]
		
		# next, remove the header by keeping on the 4th element.
		seqIDonly <- splittabs[4]
		
		# Finally, split on white space
		seqIDonly <- strsplit(seqIDonly, "[[:space:]]+")[[1]]
		
		# For each element in seqIDonly, split on the underscore delimiter
		splitseqnames <- strsplit(seqIDonly, "_", fixed=TRUE)	
		
		# make the verbose vector
		verbosesamplenamesi <- sapply(splitseqnames, function(i){i[1]})
		
		# sum the reads from each sample with tapply
		OTUi <- tapply(verbosesamplenamesi, factor(verbosesamplenamesi), length)
		
		# store results of this OTU in abundance matrix
		otumat[i, names(OTUi)] <- OTUi
	}
	
	# Return the abundance table.
	return( otu_table(otumat, taxa_are_rows=TRUE) )
}
################################################################################
#' Import new RDP OTU-table format
#' 
#' Recently updated tools on RDP Pyro site make it easier to import Pyrosequencing output 
#' into R. The modified tool ``Cluster To R Formatter'' can take a cluster file 
#' (generated from RDP Clustering tools) to create a community data matrix file
#' for distance cutoff range you are interested in. The resulting output file 
#' is a tab-delimited file containing the number of sequences for each sample 
#' for each OTU. The OTU header naming convention is \code{"OTU_"} followed by the OTU
#' number in the cluster file. It pads ``0''s to make the OTU header easy to sort.
#' The OTU numbers are not necessarily in order.
#'
#' @usage import_RDP_otu(otufile)
#'
#' @param otufile (Optional). 
#'  A character string indicating the file location of the OTU file, 
#'  produced/exported according to the instructions above.
#'
#' @return A \code{\link{otu_table-class}} object.
#'
#' @seealso
#' An alternative ``cluster'' file importer for RDP results:
#' \code{\link{import_RDP_cluster}}
#'
#' The main RDP-pyrosequencing website
#' \url{http://pyro.cme.msu.edu/index.jsp}
#'
#' @export
#' @examples
#' otufile <- system.file("extdata", "rformat_dist_0.03.txt.gz", package="phyloseq")
#' ### the gzipped file is automatically recognized, and read using R-connections
#' ex_otu  <- import_RDP_otu(otufile)
#' class(ex_otu)
#' ntaxa(ex_otu)
#' nsamples(ex_otu)
#' sample_sums(ex_otu)
#' head(t(ex_otu))
import_RDP_otu <- function(otufile){
	otumat <- read.table(otufile, TRUE, sep="\t", row.names=1)	
	return(otu_table(otumat, FALSE))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#' Imports a tab-delimited version of the pyrotagger output file.
#'
#' PyroTagger is a web-server that takes raw, barcoded 16S rRNA amplicon sequences
#' and returns an excel spreadsheet (\code{".xls"}) with both abundance and 
#' taxonomy data. It also includes some confidence information related to the 
#' taxonomic assignment.
#'
#' PyroTagger is created and maintained by the Joint Genome Institute 
#' at \code{"http://pyrotagger.jgi-psf.org/"}
#' 
#' The typical output form PyroTagger is a spreadsheet format \code{".xls"}, which poses 
#' additional import challenges. However, virtually all spreadsheet applications 
#' support the \code{".xls"} format, and can further export this file in a 
#' tab-delimited format. It is recommended that you convert the xls-file without
#' any modification (as tempting as it might be once you have loaded it) into a
#' tab-delimited text file. Deselect any options to encapsulate fields in quotes,
#' as extra quotes around each cell's contents might cause problems during 
#' file processing. These quotes will also inflate the file-size, so leave them out
#' as much as possible, while also resisting any temptation to modify the xls-file
#' ``by hand''. 
#'
#' A highly-functional and free spreadsheet application can be obtained as part
#' of the cross-platform \code{OpenOffice} suite. It works for the above 
#' required conversion. Go to \code{"http://www.openoffice.org/"}. 
#'
#' It is regrettable that this importer does not take the xls-file directly
#' as input. However, because of the moving-target nature of spreadsheet
#' file formats, there is limited support for direct import of these formats into
#' \code{R}. Rather than add to the dependency requirements of emph{phyloseq}
#' and the relative support of these xls-support packages, it seems more efficient
#' to choose an arbitrary delimited text format, and focus on the data 
#' structure in the PyroTagger output. This will be easier to support in the
#' long-run.
#' 
#' @usage import_pyrotagger_tab(pyrotagger_tab_file, 
#'	strict_taxonomy=FALSE, keep_potential_chimeras=FALSE)
#'
#' @param pyrotagger_tab_file (Required). A character string. The name of the tab-delimited
#'  pyrotagger output table.
#'
#' @param strict_taxonomy (Optional). Logical. Default \code{FALSE}. Should the taxonomyTable
#'  component be limited to just taxonomic data? Default includes all fields from 
#'  the pyrotagger file.
#'
#' @param keep_potential_chimeras (Optional). Logical. Default \code{FALSE}. The 
#'  pyrotagger output also includes OTUs that are tagged by pyrotagger as likely
#'  chimeras. These putative chimeric OTUs can be retained if set to \code{TRUE}.
#'  The putative chimeras are excluded by default.
#' 
#' @return An \code{otuTax} object containing both the otu_table and TaxonomyTable data
#'  components, parsed from the pyrotagger output.
#'
#' @export
#'
#' @references \url{http://pyrotagger.jgi-psf.org/}
#'
#' @examples
#' ## New_otuTaxObject <- import_pyrotagger_tab(pyrotagger_tab_file)
import_pyrotagger_tab <- function(pyrotagger_tab_file, 
	strict_taxonomy=FALSE, keep_potential_chimeras=FALSE){

	x <- readLines(pyrotagger_tab_file, warn=FALSE)
	# Get the header
	pyro_header <- strsplit(x[1], "\t", TRUE)[[1]]
	# Pop the first (header) line from the list.
	x <- x[-1]
	
	########################################
	### There are "Potential chimeras" 
	### listed in the typical output, separated by 2 completely blank lines
	### after the last confidently-good OTU.
	########################################
	chimera_line <- grep("Potential chimeras", x, fixed=TRUE)
	if( keep_potential_chimeras ){
		# Pop just the blank lines that delimit the chimeras
		# at the bottom of the table
		x <- x[-((chimera_line-2):chimera_line)]
	} else {
		x <- x[-((chimera_line-2):length(x))]
	}
	
	########################################
	# The tab-split character list, z
	########################################
	z <- strsplit(x, "\t", TRUE)
	names(z) <- sapply(z, function(z){z[1]})
	
	# The table switches from abundance to taxonomy at the "% Identity" column
	taxonomy_table_column_index <- which( pyro_header == "% identity" )
	
	########################################
	# Initialize the two matrices
	# (otu_table and taxonomyTable)
	########################################
	### Initialize abundance matrix, a
	a <- matrix(0, nrow=length(x), ncol=(taxonomy_table_column_index-2))
	colnames(a) <- pyro_header[2:(taxonomy_table_column_index-1)]
	rownames(a) <- names(z)
	
	###### Initialize the raw pyrotagger taxonomy matrix, w
	ntax_tablecols <- (max(sapply(z, length)) - taxonomy_table_column_index + 1)
	w <- matrix("", nrow=length(x), ncol=ntax_tablecols)
	rownames(w) <- names(z)
	colnamesw <- pyro_header[-(1:(taxonomy_table_column_index-1))]
	colnamesw <- colnamesw[1:which(colnamesw=="Taxonomy")]
	colnamesw <- c(colnamesw, paste("col", (which(colnamesw=="Taxonomy")+1):ntax_tablecols, sep="") )
	colnames(w) <- colnamesw
	
	# Rename the taxonomy columns
	biotaxonomy <- c("Domain", "Phylum", "Class", "Order",
			 	"Family", "Genus", "Species", "Strain")
	colnames(w)[which(colnames(w)=="Taxonomy"):length(colnames(w))][1:length(biotaxonomy)] <- biotaxonomy
	
	# Loop through each line and add to appropriate matrix.
	for( i in rownames(a) ){
		# i <- rownames(a)[[1]]
		# cut out just the abundance part, and convert to integer
		y <- as.integer(z[[i]][2:(taxonomy_table_column_index-1)])
		y[is.na(y)] <- 0
		a[i, ] <- y
		
		# Taxonomy data is jagged
		taxi <- z[[i]][-(1:(taxonomy_table_column_index-1))]
		w[i, 1:length(taxi)] <- taxi	
	}
	
	# Create the component objects
	OTU <- otu_table(a, taxa_are_rows=TRUE)
	if( strict_taxonomy ){
		TAX <- tax_table[, biotaxonomy]
	} else {
		TAX <- tax_table(w)
	}
	
	return( phyloseq(OTU, TAX) )

}
################################################################################
################################################################################
#' Show cutoff values available in a mothur file. 
#'
#' This is a helper function to report back to the user the different cutoff
#' values available in a given mothur file -- 
#' for instance, a list or shared file.
#'
#' @param mothur_list_file The file name and/or location as produced by \emph{mothur}.
#'
#' @return A character vector of the different cutoff values contained in the file.
#'  For a given set of arguments to the \code{cluster()} command from within
#'  \emph{mothur}, a number of OTU-clustering results are returned in the same
#'  file. The exact cutoff values used by \emph{mothur} can vary depending
#'  on the input data/parameters. This simple function returns the cutoffs that were actually
#'  included in the \emph{mothur} output. This an important extra step prior to
#'  importing data with the \code{\link{import_mothur}} function.
#' 
#' @export
#'
#' @seealso \code{\link{import_mothur}}
#'  
show_mothur_cutoffs <- function(mothur_list_file){
  unique(scan(mothur_list_file, "character", comment.char="\t", quiet=TRUE))
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
#'  file using the \code{\link{show_mothur_cutoffs}} function.
#'
#' @return A list, where each element is a character vector of 1 or more 
#'  sequence identifiers, indicating how each sequence from the original data
#'  is clustered into OTUs by \emph{mothur}. Note that in some cases this is highly
#'  dependent on the choice for \code{cutoff}. 
#'
#' @seealso \code{\link{show_mothur_cutoffs}}, \code{\link{import_mothur}}
#' @keywords internal
#'  
import_mothur_otulist <- function(mothur_list_file, cutoff=NULL){
  # mothur_list_file = system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
  # cutoff = 0.04
  cutoffs = show_mothur_cutoffs(mothur_list_file)
  cutoff = select_mothur_cutoff(cutoff, cutoffs)
  # Read only the line corresponding to that cutoff  
  inputline = which(cutoffs == cutoff)
  rawlines = scan(mothur_list_file, "character", sep="\t", skip=(inputline-1), nlines=1, na.strings="", quiet=TRUE)
  rawlines = rawlines[!is.na(rawlines)]
  # The first two elements are the cutoff and the number of OTUs. skip, and read to first comma for OTUnames
  OTUnames = scan(text=rawlines, what="character", comment.char=",", quiet=TRUE)[3:as.integer(rawlines[2])]
  # split each element on commas
  OTUs <- strsplit(rawlines[3:as.integer(rawlines[2])], ",", fixed=TRUE)
  # Name each OTU (currently as the first seq name in each cluster), and return the list
  names(OTUs) <- OTUnames
  # return as-is
  return(OTUs)
}
################################################################################
# Need to select a cutoff if none was provided by user. 
# Take the largest non-"unique" cutoff possible,
# if "unique" is the only cutoff included in the list file, use that.
# Multiple cutoffs are provided in both `.shared` and `.list` files.
# This function consolidates the heuristic for selecting/checking a specified cutoff.
#' @keywords internal
select_mothur_cutoff = function(cutoff, cutoffs){
  if( is.null(cutoff) ){
    # cutoff was NULL, need to select one.
    if( length(cutoffs) > 1 ){
      # Select the largest value, avoiding the "unique" option.
      selectCutoffs <- as(cutoffs[cutoffs != "unique"], "numeric")
      cutoff <- as.character(max(selectCutoffs))
    } else {
      # There is only one cutoff value, so use it.
      # Don't have to specify a cutoff, in this case
      cutoff <- cutoffs
    }
  } else {
    # Provided by user, non-null. Coerce to character for indexing
    cutoff <- as.character(cutoff)
    # Check that it is in set of available cutoffs.
    if( !cutoff %in% cutoffs ){
      stop("The cutoff value you provided is not among those available. Try show_mothur_cutoffs()")
    }
  }
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
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
#'
import_mothur_groups <- function(mothur_group_file){
	read.table(mothur_group_file, sep="\t", as.is=TRUE, stringsAsFactors=FALSE, colClasses="character", row.names=1)
}
################################################################################
#' Import mothur list and group files and return an otu_table
#'
#' @usage import_mothur_otu_table(mothur_list_file, mothur_group_file, cutoff=NULL)
#'
#' @param mothur_list_file The list file name and/or location as produced by \emph{mothur}.
#'
#' @param mothur_group_file The name/location of the group file produced 
#'  by \emph{mothur}'s \code{make.group()} function. It contains information
#'  about the sample source of individual sequences, necessary for creating a
#'  species/taxa abundance table (\code{otu_table}). See
#'  \code{http://www.mothur.org/wiki/Make.group}
#'
#' @param cutoff A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur} (and specified
#'  by the \code{mothur_list_file} argument). 
#'  The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_cutoffs}} function.
#'
#' @return An \code{\link{otu_table}} object.
#'
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
#' @importFrom plyr ldply
#' @importFrom plyr ddply
import_mothur_otu_table <- function(mothur_list_file, mothur_group_file, cutoff=NULL){
	otulist       <- import_mothur_otulist(mothur_list_file, cutoff)
	mothur_groups <- import_mothur_groups(mothur_group_file)
	# Initialize abundance matrix with zeros for sparse assignment
  samplenames = unique(mothur_groups[, 1])
	mothur_otu_table <- matrix(0, nrow=length(otulist), ncol=length(samplenames))
	colnames(mothur_otu_table) <- samplenames
	rownames(mothur_otu_table) <- names(otulist)

	# Write a sparse versino of the abundance table
	df = ldply(otulist, function(x){data.frame(read=x, stringsAsFactors=FALSE)})
	colnames(df)[1] <- "OTU"
	df = data.frame(df, sample=mothur_groups[df[, "read"], 1], stringsAsFactors=FALSE)
	adf = ddply(df, c("OTU", "sample"), function(x){
	  # x = subset(df, OTU=="59_3_17" & sample=="C")
	  data.frame(x[1, c("OTU", "sample"), drop=FALSE], abundance=nrow(x))
	})
	
	# Vectorized for speed using matrix indexing.
	# See help("Extract") for details about matrix indexing. Diff than 2-vec index.
	mothur_otu_table[as(adf[, c("OTU", "sample")], "matrix")] <- adf[, "abundance"] 
	
	# Finally, return the otu_table as a phyloseq otu_table object.
	return(otu_table(mothur_otu_table, taxa_are_rows=TRUE))
}
################################################################################
#' Import mothur shared file and return an otu_table
#'
#' @param mothur_shared_file (Required). A 
#' \href{http://www.mothur.org/wiki/Shared_file}{shared file}
#' produced by \emph{mothur}.
#'
#' @return An \code{\link{otu_table}} object.
#'
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
import_mothur_shared = function(mothur_shared_file, cutoff=NULL){
  #mothur_shared_file = "~/github/phyloseq/inst/extdata/esophagus.fn.shared.gz"
  # Check that cutoff is in cutoffs, or select a cutoff if none given.
  cutoffs = show_mothur_cutoffs(mothur_shared_file)
  cutoffs = cutoffs[!cutoffs %in% "label"]
  cutoff = select_mothur_cutoff(cutoff, cutoffs)
  x = readLines(mothur_shared_file)
  rawtab = read.table(text=x[grep(paste0("^", cutoff), x)], header=FALSE, row.names=2, stringsAsFactors=FALSE)[, -(1:2)]
  colnames(rawtab) <- strsplit(x[1], "\t")[[1]][4:(ncol(rawtab)+3)]
  return(otu_table(t(as.matrix(rawtab)), taxa_are_rows=TRUE))
}
################################################################################
#' Import mothur constaxonomy file and return a taxonomyTable
#'
#' @param mothur_constaxonomy_file (Required). A 
#'  \href{http://www.mothur.org/wiki/Constaxonomy_file}{consensus taxonomy file} 
#'  produced by \emph{mothur}.
#'  
#' @param parseFunction (Optional). A specific function used for parsing the taxonomy string.
#'  See \code{\link{parse_taxonomy_default}} for an example. If the default is
#'  used, this function expects a semi-colon delimited taxonomy string, with
#'  no additional rank specifier. A common taxonomic database is GreenGenes,
#'  and for recent versions its taxonomy includes a prefix, which is best cleaved
#'  and used to precisely label the ranks (\code{\link{parse_taxonomy_greengenes}}).
#'
#' @return An \code{\link{taxonomyTable-class}} object.
#'
#' @seealso \code{\link{import_mothur}}
#' 
#' \code{\link{tax_table}}
#' 
#' \code{\link{phyloseq}}
#' 
#' @keywords internal
import_mothur_constaxonomy = function(mothur_constaxonomy_file, parseFunction=parse_taxonomy_default){
  read.table(mothur_constaxonomy_file)
  rawtab = read.table(mothur_constaxonomy_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)[, "Taxonomy", drop=FALSE]
  if( identical(parseFunction, parse_taxonomy_default) ){
    # Proceed with default parsing stuff.
    # Remove the confidence strings inside the parentheses, if present
    rawtab[, "Taxonomy"] = gsub("\\([[:digit:]]+\\)", "", rawtab[, "Taxonomy"])
    # Remove the quotation marks, if present
    rawtab[, "Taxonomy"] = gsub("\"", "", rawtab[, "Taxonomy"])
    # Remove trailing semicolon
    rawtab[, "Taxonomy"] = gsub(";$", "", rawtab[, "Taxonomy"])
    # Split on semicolon
    taxlist = strsplit(rawtab[, "Taxonomy"], ";", fixed=TRUE)
    taxlist = lapply(taxlist, parseFunction)
  } else {
    taxlist = lapply(rawtab[, "Taxonomy"], parseFunction)
  }
  names(taxlist) <- rownames(rawtab)
  return(build_tax_table(taxlist))
}
################################################################################
#' General function for importing mothur data files into phyloseq.
#'
#' Technically all parameters are optional,
#' but if you don't provide any file connections, then nothing will be returned.
#' While the \code{list} and \code{group} files are the first two arguments
#' for legacy-compatibility reasons, we don't recommend that you use these
#' file types with modern (large) datasets. They are comically inefficient, as
#' they store the name of every sequencing read in both files. The \emph{mothur}
#' package provides conversions utilities to create other more-efficient formats,
#' which we recommend, like 
#' the \href{http://www.mothur.org/wiki/Shared_file}{shared file} for an OTU table.
#' Alternatively, mothur also provides a utility to create a biom-format file
#' that is independent of OTU clustering platform. Biom-format files 
#' should be imported not with this function, but with \code{\link{import_biom}}.
#' The resulting objects after import should be \code{\link{identical}} in R.
#'
#' @param mothur_list_file (Optional). The list file name / location produced by \emph{mothur}.
#'
#' @param mothur_group_file (Optional). The name/location of the group file produced 
#'  by \emph{mothur}'s \code{make.group()} function. It contains information
#'  about the sample source of individual sequences, necessary for creating a
#'  species/taxa abundance table (\code{otu_table}). See
#'  \code{http://www.mothur.org/wiki/Make.group} 
#'
#' @param mothur_tree_file (Optional). 
#'  A tree file, presumably produced by \emph{mothur},
#'  and readable by \code{\link{read_tree}}.
#'  The file probably has extension \code{".tree"}.
#'
#' @param cutoff (Optional). A character string indicating the cutoff value, (or \code{"unique"}), 
#'  that matches one of the cutoff-values used to produce the OTU clustering 
#'  results contained within the list-file created by \emph{mothur} (and specified
#'  by the \code{mothur_list_file} argument). The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_cutoffs}} function.
#'
#' @param mothur_shared_file (Optional). A 
#' \href{http://www.mothur.org/wiki/Shared_file}{shared file}
#' produced by \emph{mothur}.
#' 
#' @param mothur_constaxonomy_file (Optional). A 
#'  \href{http://www.mothur.org/wiki/Constaxonomy_file}{consensus taxonomy file} 
#'  produced by \emph{mothur}.
#'  
#' @param parseFunction (Optional). A specific function used for parsing the taxonomy string.
#'  See \code{\link{parse_taxonomy_default}} for an example. If the default is
#'  used, this function expects a semi-colon delimited taxonomy string, with
#'  no additional rank specifier. A common taxonomic database is GreenGenes,
#'  and in recent versions its taxonomy entries include a prefix, which is best cleaved
#'  and used to precisely label the ranks (\code{\link{parse_taxonomy_greengenes}}).
#'
#' @return The object class depends on the provided arguments. 
#'  A phyloseq object is returned if enough data types are provided.
#'  If only one data component can be created from the data, it is returned.
#' 
#'  FASTER (recommended for larger data sizes):
#'  
#'  If only a \code{mothur_constaxonomy_file} is provided, 
#'  then a  \code{\link{taxonomyTable-class}} object is returned.
#'  
#'  If only a \code{mothur_shared_file} is provided, 
#'  then an \code{\link{otu_table}} object is returned.
#'  
#'  SLOWER (but fine for small file sizes):
#'  
#'  The list and group file formats are extremely inefficient for large datasets,
#'  and they are not recommended. The mothur software provides tools for 
#'  converting to other file formats, such as a so-called ``shared'' file.
#'  You should provide a shared file, or group/list files, but not
#'  both at the same time.
#'  If only a list and group file are provided, 
#'  then an \code{otu_table} object is returned. 
#'  Similarly, if only a list and tree file are provided, 
#'  then only a tree is returned (\code{\link[ape]{phylo}}-class).
#'
#' @references \url{http://www.mothur.org/wiki/Main_Page}
#'
#' Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent,
#' community-supported software for describing and comparing microbial communities.
#' Appl Environ Microbiol, 2009. 75(23):7537-41.
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
#' # show_mothur_cutoffs(mothur_list_file)
#' # test1 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file)
#' # test2 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file, cutoff="0.02")
#' # # Returns just a tree
#' # import_mothur(mothur_list_file, mothur_tree_file=mothur_tree_file)
#' # # Returns just an otu_table
#' # import_mothur(mothur_list_file, mothur_group_file=mothur_group_file)
#' # # Returns an error
#' # import_mothur(mothur_list_file)
#' # # Should return an "OMG, you must provide the list file" error
#' # import_mothur()
import_mothur <- function(mothur_list_file=NULL,  mothur_group_file=NULL,
  mothur_tree_file=NULL, cutoff=NULL,
  mothur_shared_file=NULL, mothur_constaxonomy_file=NULL, parseFunction=parse_taxonomy_default){

  pslist = vector("list")
  
	if( !is.null(mothur_group_file) & !is.null(mothur_list_file) ){
	  # If list & group files provided, you can make an OTU table.
		groupOTU = import_mothur_otu_table(mothur_list_file, mothur_group_file, cutoff)
		pslist = c(pslist, list(groupOTU))
	} 
	
	if( !is.null(mothur_tree_file) ){
		tree <- read_tree(mothur_tree_file)
		pslist = c(pslist, list(tree))
	}
  
  if( !is.null(mothur_shared_file) ){
    OTUshared <- import_mothur_shared(mothur_shared_file)
    pslist = c(pslist, list(OTUshared))
  }
  
  if( !is.null(mothur_constaxonomy_file) ){
    tax <- import_mothur_constaxonomy(mothur_constaxonomy_file, parseFunction)
    pslist = c(pslist, list(tax))
  }  
  
  return(do.call("phyloseq", pslist))
}
################################################################################
#' Import mothur-formatted distance file
#'
#' The mothur application will produce a file containing the pairwise distances
#' between all sequences in a dataset. This distance matrix can be the basis for
#' OTU cluster designations. R also has many built-in or off-the-shelf tools for
#' dealing with distance matrices. 
#' 
#' @usage import_mothur_dist(mothur_dist_file)
#' 
#' @param mothur_dist_file Required. The distance file name / location produced by \emph{mothur}.
#'
#' @return A distance matrix object describing all sequences in a dataset.
#'
#' @export
#'
#' @seealso \code{\link{import_mothur}}
#'
#' @examples
#' # # Take a look at the dataset shown here as an example:
#' # # "http://www.mothur.org/wiki/Esophageal_community_analysis"
#' # # find the file ending with extension ".dist", download to your system
#' # # The location of your file may vary
#' # mothur_dist_file <- "~/Downloads/mothur/Esophagus/esophagus.dist"
#' # myNewDistObject  <- import_mothur_dist(mothur_dist_file)
import_mothur_dist <- function(mothur_dist_file){
	# Read the raw distance matrix file produced by mothur:
	raw_dist_lines <- readLines(mothur_dist_file)
	
	# split each line on white space, and begin modifying into dist-matrix format
	dist_char <- strsplit(raw_dist_lines, "[[:space:]]+")
	dist_char <- dist_char[-1]
	# add name to each list element
	names(dist_char) <- sapply(dist_char, function(i){i[1]})
	# pop out the names from each vector
	dist_char <- sapply(dist_char, function(i){i[-1]})
	# convert to numeric vectors
	dist_char <- sapply(dist_char, function(i){as(i, "numeric")})
	
	# Initialize and fill the matrix
	distm <- matrix(0, nrow=length(dist_char), ncol=length(dist_char))
	rownames(distm) <- names(dist_char); colnames(distm) <- names(dist_char)
	for( i in names(dist_char)[-1] ){
		distm[i, 1:length(dist_char[[i]])] <- dist_char[[i]]
	}
	diag(distm) <- 1
	distd <- as.dist(distm)
	return(distd)
}
################################################################################
################################################################################
#' Export a distance object as \code{.names} and \code{.dist} files for mothur
#'
#' The purpose of this function is to allow a user to easily export a distance object
#' as a pair of files that can be immediately imported by mothur for OTU clustering
#' and related analysis. A distance object can be created in \code{R} in a number of
#' ways, including via cataloguing the cophentic distances of a tree object.
#'
#' @usage export_mothur_dist(x, out=NULL, makeTrivialNamesFile=NULL)
#'
#' @param x (Required). A \code{"dist"} object, or a symmetric matrix.
#'
#' @param out (Optional). The desired output filename for the \code{.dist} file, OR
#'  left \code{NULL}, the default, in which case the mothur-formated distance table
#'  is returned to \code{R} standard out.
#'
#' @param makeTrivialNamesFile (Optional). Default \code{NULL}. The desired name of the \code{.names} file.
#'  If left \code{NULL}, the file name will be a modified version of the \code{out} argument.
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
#' @examples #
#' data(esophagus) 
#' myDistObject <- as.dist(ape::cophenetic.phylo(phy_tree(esophagus)))
#' export_mothur_dist(myDistObject)
export_mothur_dist <- function(x, out=NULL, makeTrivialNamesFile=NULL){
	if( class(x)== "matrix" ){ x <- as.dist(x) }
	if( class(x)!= "dist" ){ stop("x must be a dist object, or symm matrix") }

	# While x is a dist-object, get the length of unique pairs
	# to initialize the dist table.
	distdf <- matrix("", nrow=length(x), ncol=3)

	# Now convert x to matrix for looping, indexing.
	x <- as(x, "matrix")
	colnames(distdf) <- c("i", "j", "d") 
	# distdf row counter
	z <- 1
	
	# The big loop.
	for( i in 2:nrow(x) ){ # i <- 2
		thisvec <- x[i, 1:(i-1)]
		for( j in 1:length(thisvec) ){ # j <- 1
			distdf[z, "i"] <- rownames(x)[i]
			distdf[z, "j"] <- colnames(x)[j]
			distdf[z, "d"] <- thisvec[j]
			z <- z + 1
		}
	}	

	# mothur requires a .names file, in case you removed identical sequences
	# from within mothur and need to keep track and add them back.
	if( !is.null(makeTrivialNamesFile) ){
		namestab <- matrix(rownames(x), nrow=length(rownames(x)), ncol=2)
		write.table(namestab, file=makeTrivialNamesFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
	
	# If is.null(out)==TRUE, then return two-column table.
	# If it's a character, write.table-it 
	if( is.null(out) ){
		return(distdf)
	} else {
		write.table(distdf, file=out, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
}
################################################################################
#' Export environment (ENV) file for UniFrac Server.
#'
#' Creates the environment table that is needed for the original UniFrac
#' algorithm. Useful for cross-checking, or if want to use UniFrac server.
#' Optionally the ENV-formatted table can be returned to the \code{R}
#' workspace, and the tree component can be exported as Nexus format
#' (Recommended). 
#'
#' @param physeq (Required). Experiment-level (\code{\link{phyloseq-class}}) object.
#'  Ideally this also contains the phylogenetic tree, which is also exported by default.
#'
#' @param file (Optional). The file path for export. If not-provided, the 
#'  expectation is that you will want to set \code{return} to \code{TRUE},
#'  and manipulate the ENV table on your own. Default is \code{""}, skipping 
#'  the ENV file from being written to a file.
#'
#' @param writeTree (Optional). Write the phylogenetic tree as well as the 
#'  the ENV table. Default is \code{TRUE}.
#'
#' @param return (Optional). Should the ENV table be returned to the R workspace?
#'  Default is \code{FALSE}.
#'
#' @importFrom ape write.nexus
#' @export
#' 
#' @examples
#' # # Load example data
#' # data(esophagus)
#' # export_env_file(esophagus, "~/Desktop/esophagus.txt")
export_env_file <- function(physeq, file="", writeTree=TRUE, return=FALSE){
	# data(esophagus)
	# physeq <- esophagus

	# Create otu_table matrix and force orientation
	OTU     <- as(otu_table(physeq), "matrix")
	if( !taxa_are_rows(physeq) ){ OTU <- t(OTU) }
	
	# initialize sequence/sample names
	seqs    <- taxa_names(physeq)
	samples <- sample_names(physeq)	
	
	# initialize output table as matrix
	ENV <- matrix("", nrow=sum(OTU >= 1), ncol=3)
	
	# i counts the row of the output , ENV
	i=1
	while( i < nrow(ENV) ){
		for( j in seqs){ 
			for( k in which(OTU[j, ]>0) ){
					ENV[i, 1] <- j
					ENV[i, 2] <- samples[k]
					ENV[i, 3] <- OTU[j, k]
					i <- i + 1
			}
		}
	}
	# If a file path is provided, write the table to that file
	if(file != ""){
		write.table(ENV, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
	
	# If needed, also write the associated tree-file. 
	if( writeTree ){
		fileTree <- paste(file, ".nex", sep="")
		write.nexus(phy_tree(physeq), file=fileTree, original.data=FALSE)
	}

	# If return argument is TRUE, return the environment table
	if(return){ return(ENV) }
}
################################################################################
# UniFrac ENV files have the form:
# 
# SEQ1        ENV1        1
# SEQ1        ENV2        2
# SEQ2        ENV1        15
# SEQ3        ENV1        2
# SEQ4        ENV2        8
# SEQ5        ENV1        4
# http://128.138.212.43/unifrac/help.psp#env_file
################################################################################
#' Import phyloseq data from biom-format file
#'
#' New versions of QIIME produce a more-comprehensive and formally-defined
#' JSON file format, called biom file format:
#'
#' ``The biom file format (canonically pronounced `biome') is designed to be a 
#' general-use format for representing counts of observations in one or
#' more biological samples. BIOM is a recognized standard for the Earth Microbiome
#' Project and is a Genomics Standards Consortium candidate project.''
#'
#' \url{http://biom-format.org/} 
#'
#' @usage import_biom(BIOMfilename, 
#'  treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
#'  parseFunction=parse_taxonomy_default, parallel=FALSE, version=1.0, ...)
#'
#' @param BIOMfilename (Required). A character string indicating the 
#'  file location of the BIOM formatted file. This is a JSON formatted file,
#'  specific to biological datasets, as described in 
#'  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{the biom-format home page}.
#'  In principle, this file should include you OTU abundance data (OTU table),
#'  your taxonomic classification data (taxonomy table), as well as your
#'  sample data, for instance what might be in your ``sample map'' in QIIME.
#'  A phylogenetic tree is not yet supported by biom-format, and so is a
#'  separate argument here. If, for some reason, your biom-format file is
#'  missing one of these mentioned data types but you have it in a separate file,
#'  you can first import the data that is in the biom file using this function,
#'  \code{import_biom}, and then ``merge'' the remaining data after you have
#'  imported with other tools using the relatively general-purpose data
#'  merging function called \code{\link{merge_phyloseq}}.
#'
#' @param treefilename (Optional). Default value is \code{NULL}. 
#'  A file representing a phylogenetic tree
#'  or a \code{\link{phylo}} object.
#'  Files can be NEXUS or Newick format.
#'  See \code{\link{read_tree}} for more details.
#'  Also, if using a recent release of the GreenGenes database tree,
#'  try the \code{\link{read_tree_greengenes}} function --
#'  this should solve some issues specific to importing that tree. 
#'  If provided, the tree should have the same OTUs/tip-labels
#'  as the OTUs in the other files. 
#'  Any taxa or samples missing in one of the files is removed from all. 
#'  As an example from the QIIME pipeline,
#'  this tree would be a tree of the representative 16S rRNA sequences from each OTU
#'  cluster, with the number of leaves/tips equal to the number of taxa/species/OTUs,
#'  or the complete reference database tree that contains the OTU identifiers
#'  of every OTU in your abundance table.
#'  Note that this argument can be a tree object (\code{\link[ape]{phylo}}-class)
#'  for cases where the tree has been --- or needs to be --- imported separately,
#'  as in the case of the GreenGenes tree mentioned earlier (code{\link{read_tree_greengenes}}).
#'
#' @param refseqfilename (Optional). Default \code{NULL}.
#'  The file path of the biological sequence file that contains at a minimum
#'  a sequence for each OTU in the dataset.
#'  Alternatively, you may provide an already-imported
#'  \code{\link[Biostrings]{XStringSet}} object that satisfies this condition.
#'  In either case, the \code{\link{names}} of each OTU need to match exactly the
#'  \code{\link{taxa_names}} of the other components of your data.
#'  If this is not the case, for example if the data file is a FASTA format but
#'  contains additional information after the OTU name in each sequence header,
#'  then some additional parsing is necessary,
#'  which you can either perform separately before calling this function, 
#'  or describe explicitly in a custom function provided in the (next) argument,
#'  \code{refseqFunction}.
#'  Note that the \code{\link[Biostrings]{XStringSet}} class can represent any 
#'  arbitrary sequence, including user-defined subclasses, but is most-often
#'  used to represent RNA, DNA, or amino acid sequences. 
#'  The only constraint is that this special list of sequences
#'  has exactly one named element for each OTU in the dataset.
#' 
#' @param refseqFunction (Optional). 
#'  Default is \code{\link[Biostrings]{readDNAStringSet}}, 
#'  which expects to read a fasta-formatted DNA sequence file.
#'  If your reference sequences for each OTU are amino acid, RNA, or something else,
#'  then you will need to specify a different function here.
#'  This is the function used to read the file connection provided as the
#'  the previous argument, \code{refseqfilename}.
#'  This argument is ignored if \code{refseqfilename} is already a
#'  \code{\link[Biostrings]{XStringSet}} class. 
#' 
#' @param refseqArgs (Optional). 
#'  Default \code{NULL}.
#'  Additional arguments to \code{refseqFunction}.
#'  See \code{\link[Biostrings]{XStringSet-io}} for details about
#'  additional arguments to the standard read functions in the Biostrings package.
#'
#' @param parseFunction (Optional). A function. It must be a function that
#'  takes as its first argument a character vector of taxonomic rank labels
#'  for a single OTU
#'  and parses and names each element
#'  (an optionally removes unwanted elements).
#'  Further details and examples of acceptable functions are provided
#'  in the documentation for \code{\link{parse_taxonomy_default}}. 
#'  There are many variations on taxonomic nomenclature, and naming 
#'  conventions used to store that information in various taxonomic 
#'  databases and phylogenetic assignment algorithms. A popular database,
#'  \url{http://greengenes.lbl.gov/cgi-bin/nph-index.cgi}{greengenes},
#'  has its own custom parsing function provided in the phyloseq package,
#'  \code{\link{parse_taxonomy_greengenes}},
#'  and more can be contributed or posted as code snippets as needed.
#'  They can be custom-defined by a user immediately prior to the the call to
#'  \code{\link{import_biom}}, and this is a suggested first step to take
#'  when trouble-shooting taxonomy-related errors during file import.
#'
#' @param parallel (Optional). Logical. Wrapper option for \code{.parallel}
#'  parameter in \code{plyr-package} functions. If \code{TRUE}, apply 
#'  parsing functions in parallel, using parallel backend provided by
#'  \code{\link{foreach}} and its supporting backend packages. One caveat,
#'  plyr-parallelization currently works most-cleanly with \code{multicore}-like
#'  backends (Mac OS X, Unix?), and may throw warnings for SNOW-like backends.
#'  See the example below for code invoking multicore-style backend within
#'  the \code{doParallel} package.
#'
#'  Finally, for many datasets a parallel import should not be necessary
#'  because a serial import will be just as fast and the import is often only
#'  performed one time; after which the data should be saved as an RData file
#'  using the \code{\link{save}} function.
#' 
#' @param version (Optional). Numeric. The expected version number of the file.
#'  As the BIOM format evolves, version-specific importers may be available
#'  by adjusting the version value. Default is \code{1.0}. 
#'  Not yet implemented. Parsing of the biom-format is done mostly
#'  by the biom package now available in CRAN.
#'
#' @param ... Additional parameters passed on to \code{\link{read_tree}}.
#'
#' @return A \code{\link{phyloseq-class}} object.
#'
#' @seealso
#' \code{\link{import}}
#'
#' \code{\link{import_qiime}}
#'
#' \code{\link{read_tree}}
#' 
#' \code{\link{read_tree_greengenes}}
#' 
#' \code{\link[biom]{biom-package}}
#' 
#' \code{\link[biom]{read_biom}}
#'
#' \code{\link[biom]{biom_data}}
#'
#' \code{\link[biom]{sample_metadata}}
#' 
#' \code{\link[biom]{observation_metadata}}
#'
#' \code{\link[Biostrings]{XStringSet-io}}
#'
#' @references \href{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{biom-format}
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom biom read_biom
#' @importFrom biom sample_metadata
#' @importFrom biom biom_data
#' @importFrom biom observation_metadata
#' @export
#' @examples
#' # An included example of a rich dense biom file
#' rich_dense_biom  <- system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
#' import_biom(rich_dense_biom,  parseFunction=parse_taxonomy_greengenes)
#' # An included example of a sparse dense biom file
#' rich_sparse_biom <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
#' import_biom(rich_sparse_biom, parseFunction=parse_taxonomy_greengenes)
#' # # # Example code for importing large file with parallel backend
#' # library("doParallel")
#' # registerDoParallel(cores=6)
#' # import_biom("my/file/path/file.biom", parseFunction=parse_taxonomy_greengenes, parallel=TRUE)
import_biom <- function(BIOMfilename, 
	treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
	parseFunction=parse_taxonomy_default, parallel=FALSE, version=1.0, ...){

	# initialize the argument-list for phyloseq. Start empty.
	argumentlist <- list()
	
	# Read the data
	x = read_biom(biom_file=BIOMfilename)
	
	########################################
	# OTU table:
	########################################
	otutab = otu_table(as(biom_data(x), "matrix"), taxa_are_rows=TRUE)
	argumentlist <- c(argumentlist, list(otutab))
	
	########################################
	# Taxonomy Table
	########################################
	# Need to check if taxonomy information is empty (minimal BIOM file)
	if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
	  taxtab <- NULL
  } else {
    # parse once each character vector, save as a list
    taxlist = lapply(x$rows, function(i){
      parseFunction(i$metadata$taxonomy)
    })
    names(taxlist) = sapply(x$rows, function(i){i$id})
    taxtab = build_tax_table(taxlist)
	}
	argumentlist <- c(argumentlist, list(taxtab))
	
	########################################
	# Sample Data ("columns" in QIIME/BIOM)
	########################################
	# If there is no metadata (all NULL), then set sam_data <- NULL
	if( is.null(sample_metadata(x)) ){
		samdata <- NULL
	} else {
		samdata = sample_data(sample_metadata(x))
	}
	argumentlist <- c(argumentlist, list(samdata))

	########################################
	# Tree data
	########################################
	if( !is.null(treefilename) ){
		if( inherits(treefilename, "phylo") ){
			# If argument is already a tree, don't read, just assign.
			tree = treefilename
		} else {
			# NULL is silently returned if tree is not read properly.
			tree <- read_tree(treefilename, ...)		
		}
		# Add to argument list or warn
		if( is.null(tree) ){
			warning("treefilename failed import. It not included.")
		} else {
			argumentlist <- c(argumentlist, list(tree) )
		}
	}
	
	########################################
	# Reference Sequence data
	########################################	
	if( !is.null(refseqfilename) ){
		if( inherits(refseqfilename, "XStringSet") ){
			# If argument is already a XStringSet, don't read, just assign.
			refseq = refseqfilename
		} else {
			# call refseqFunction and read refseqfilename, either with or without additional args
			if( !is.null(refseqArgs) ){
				refseq = do.call("refseqFunction", c(list(refseqfilename), refseqArgs))
			} else {
				refseq = refseqFunction(refseqfilename)
			}
		}
		argumentlist <- c(argumentlist, list(refseq) )
	}
	
	########################################
	# Put together into a phyloseq object
	########################################
	return( do.call("phyloseq", argumentlist) )

}
################################################################################
# Need to export these parsing functions as examples...
################################################################################
#' Parse elements of a taxonomy vector
#'
#' These are provided as both example and default functions for
#' parsing a character vector of taxonomic rank information for a single taxa.
#' As default functions, these are intended for cases where the data adheres to
#' the naming convention used by greengenes
#' (\url{http://greengenes.lbl.gov/cgi-bin/nph-index.cgi})
#' or where the convention is unknown, respectively.
#' To work, these functions -- and any similar custom function you may want to
#' create and use -- must take as input a single character vector of taxonomic
#' ranks for a single OTU, and return a \strong{named} character vector that has
#' been modified appropriately (according to known naming conventions,
#' desired length limits, etc.
#' The length (number of elements) of the output named vector does \strong{not}
#' need to be equal to the input, which is useful for the cases where the
#' source data files have extra meaningless elements that should probably be
#' removed, like the ubiquitous 
#' ``Root'' element often found in greengenes/QIIME taxonomy labels.
#' In the case of \code{parse_taxonomy_default}, no naming convention is assumed and
#' so dummy rank names are added to the vector.  
#' More usefully if your taxonomy data is based on greengenes, the
#' \code{parse_taxonomy_greengenes} function clips the first 3 characters that 
#' identify the rank, and uses these to name the corresponding element according
#' to the appropriate taxonomic rank name used by greengenes
#' (e.g. \code{"p__"} at the beginning of an element means that element is 
#' the name of the phylum to which this OTU belongs).
#' Most importantly, the expectations for these functions described above
#' make them compatible to use during data import,
#' specifcally the \code{\link{import_biom}} function, but 
#' it is a flexible structure that will be implemented soon for all phyloseq
#' import functions that deal with taxonomy (e.g. \code{\link{import_qiime}}).
#'
#' @usage parse_taxonomy_default(char.vec)
#' @usage parse_taxonomy_greengenes(char.vec)
#' @usage parse_taxonomy_qiime(char.vec)
#'
#' @param char.vec (Required). A single character vector of taxonomic
#'  ranks for a single OTU, unprocessed (ugly).
#' 
#' @return A character vector in which each element is a different
#'  taxonomic rank of the same OTU, and each element name is the name of
#'  the rank level. For example, an element might be \code{"Firmicutes"}
#'  and named \code{"phylum"}.
#'  These parsed, named versions of the taxonomic vector should
#'  reflect embedded information, naming conventions,
#'  desired length limits, etc; or in the case of \code{\link{parse_taxonomy_default}},
#'  not modified at all and given dummy rank names to each element.
#'
#' @rdname parseTaxonomy-functions
#' @export
#'
#' @seealso
#'  \code{\link{import_biom}}
#'  \code{\link{import_qiime}}
#'
#' @examples
#'  taxvec1 = c("Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae")
#'  parse_taxonomy_default(taxvec1)
#'  parse_taxonomy_greengenes(taxvec1)
#'  taxvec2 = c("Root;k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae")
#'  parse_taxonomy_qiime(taxvec2)
parse_taxonomy_default = function(char.vec){
	# Remove any leading empty space
	char.vec = gsub("^[[:space:]]{1,}", "", char.vec)
	# Remove any trailing space
	char.vec = gsub("[[:space:]]{1,}$", "", char.vec)
	if( length(char.vec) > 0 ){
		# Add dummy element (rank) name
		names(char.vec) = paste("Rank", 1:length(char.vec), sep="")
	} else {
		warning("Empty taxonomy vector encountered.")
	}
	return(char.vec)
}
#' @rdname parseTaxonomy-functions
#' @aliases parse_taxonomy_default
#' @export
parse_taxonomy_greengenes <- function(char.vec){
	# Use default to assign names to elements in case problem with greengenes prefix
	char.vec = parse_taxonomy_default(char.vec)
	# Define the meaning of each prefix according to GreenGenes taxonomy
	Tranks = c(k="Kingdom", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species")
	# Check for prefix using regexp, warn if there were none. trim indices, ti
	ti = grep("[[:alpha:]]{1}\\_\\_", char.vec)
	if( length(ti) == 0L ){
		warning(
			"No greengenes prefixes were found. \n",
			"Consider using parse_taxonomy_default() instead if true for all OTUs. \n",
			"Dummy ranks may be included among taxonomic ranks now."
		)
		# Will want to return without further modifying char.vec
		taxvec = char.vec
	# Replace names of taxvec according to prefix, if any present...
	} else {
		# Remove prefix using sub-"" regexp, call result taxvec
		taxvec = gsub("[[:alpha:]]{1}\\_\\_", "", char.vec)
		# Define the ranks that will be replaced
		repranks = Tranks[substr(char.vec[ti], 1, 1)]
		# Replace, being sure to avoid prefixes not present in Tranks
		names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
	}
	return(taxvec)
}
#' @rdname parseTaxonomy-functions
#' @aliases parse_taxonomy_default
#' @export
parse_taxonomy_qiime <- function(char.vec){
	parse_taxonomy_greengenes(strsplit(char.vec, ";", TRUE)[[1]])
}
################################################################################
#' Build a \code{\link{tax_table}} from a named possibly-jagged list
#'
#' @param taxlist (Required). A list in which each element is a vector of 
#'  taxonomic assignments named by rank.
#'  Every element of every vector must be named by the rank it represents.
#'  Every element of the list (every vector) should correspond to a single OTU
#'  and be named for that OTU. 
#'
#' @return A \code{\link{tax_table}} (\code{\link{taxonomyTable-class}}) that 
#'  has been built from \code{taxlist}. The OTU names of this output will be
#'  the element names of \code{taxlist}, and a separate taxonomic rank
#'  (column) will be included for each unique rank found among the element names
#'  of each vector in the list. \code{NA_character_} is the default value of
#'  elements in the \code{\link{tax_table}} for which there is no corresponding
#'  information in \code{taxlist}.
#'
#' @seealso
#'  \code{\link{import_biom}}
#'  \code{\link{import_qiime}}
#'
#' @export
#'
#' @examples
#'  taxvec1 = c("Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae")
#'  parse_taxonomy_default(taxvec1)
#'  parse_taxonomy_greengenes(taxvec1)
#'  taxvec2 = c("Root;k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae")
#'  parse_taxonomy_qiime(taxvec2)
#'  taxlist1 = list(OTU1=parse_taxonomy_greengenes(taxvec1), OTU2=parse_taxonomy_qiime(taxvec2))
#'  taxlist2 = list(OTU1=parse_taxonomy_default(taxvec1), OTU2=parse_taxonomy_qiime(taxvec2))
#'  build_tax_table(taxlist1)
#'  build_tax_table(taxlist2)
build_tax_table = function(taxlist){
	# Determine column headers (rank names) of taxonomy table
	columns = unique(unlist(lapply(taxlist, names)))
	# Initialize taxonomic character matrix
	taxmat <- matrix(NA_character_, nrow=length(taxlist), ncol=length(columns))
	colnames(taxmat) = columns
	# Fill in the matrix by row.
	for( i in 1:length(taxlist) ){
		# Protect against empty taxonomy
		if( length(taxlist[[i]]) > 0 ){
			# The extra column name check solves issues with raggedness, and disorder.
			taxmat[i, names(taxlist[[i]])] <- taxlist[[i]]
		}
	}
	# Convert functionally empty elements, "", to NA
	taxmat[taxmat==""] <- NA_character_
	# Now coerce to matrix, name the rows as "id" (the taxa name), coerce to taxonomyTable
	taxmat			 <- as(taxmat, "matrix")
	rownames(taxmat) = names(taxlist)
	return( tax_table(taxmat) )
}
################################################################################
################################################################################
################################################################################
#' Download and import directly from microbio.me/qiime
#' 
#' This function is for accessing microbiome datasets from the
#' \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}
#' public repository from within R.
#' I haven't yet figured out how to list the available studies at
#' \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}
#' from within R
#' (and I don't know if they have made this possible),
#' but I have provided illustrated instructions for finding 
#' details about studies you might want to download and explore,
#' as well as the FTP address that you will need. 
#' Note that you likely need to create an account at
#' \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}
#' in order to explore the studies they have available for download.
#' Please see a detailed tutorial for further help finding the FTP address 
#' of the data that you want. Also note that this function will work
#' on the full address to any \code{.zip} zipped file with microbiome
#' data organized with the same naming conventions and file types as
#' \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}.
#' Please also see the 
#' \href{http://joey711.github.io/phyloseq/download-microbio.me.html}{microbio_me_qiime tutorial}
#' 
#' @param zipftp (Required). A character string that is the full URL
#'  path to a zipped file that follows the file naming conventions used by 
#'  \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}.
#'  Alternatively, you can simply provide the study number,
#'  as a single-length \code{\link{numeric}},
#'  and this function will complete the remainder of the ftp URL hosted at
#'  \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}.
#'  For example, instead of the full URL string, 
#'  \code{"ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_494_split_library_seqs_and_mapping.zip"}, you could simply provide \code{494} as the `zipftp` argument.
#'  Note that this is internally dispatching based on the class,
#'  so \code{"494"} will not work (just leave off the quotes for study number-only).
#'  
#'  @param ext (Optional). A \code{\link{character}} string of the expected
#'   file extension, which also indicates the compression type,
#'   if \code{zipftp} is a study number instead of the full path. 
#'   Note that this argument has no effect if \code{zipftp} is the full path,
#'   in which case the file extension is read directly from the downloaded file.
#'   
#'  @param parsef (Optional). The type of taxonomic parsing to use for the
#'   OTU taxonomic classification, in the \code{.biom} file, if present.
#'   This is passed on to \code{\link{import_biom}}, but unlike that function
#'   the default parsing function is \code{\link{parse_taxonomy_greengenes}},
#'   rather than \code{\link{parse_taxonomy_default}}, because we know
#'   ahead of time that most (or all?) of the taxonomic classifications
#'   in the \code{microbio.me/qiime} repository will be based on 
#'   GreenGenes.
#'  
#'  @param ... (Optional, for advanced users). Additional arguments passed to 
#'   \code{\link{download.file}}. This is mainly for non-standard links to
#'   resources (in this case, a zipped file) that are not being hosted by
#'   \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}.
#'   If you are using a FTP address or study number from their servers,
#'   then you shouldn't need to provide any additional arguments.
#' 
#' @return
#'  A \code{\link{phyloseq-class}} object if possible, a component if only a 
#'  component could be imported, or \code{NULL} if nothing could be imported
#'  after unzipping the file. Keep in mind there is a specific naming-convention
#'  that is expected based on the current state of the
#'  \href{http://www.microbio.me/qiime/index.psp}{microbio.me/qiime}
#'  servers. Several helpful messages are \code{\link{cat}}ted to standard out
#'  to help let you know the ongoing status of the current 
#'  download and import process.
#'
#' @seealso
#'  See \code{\link{download.file}} and \code{\link{url}}
#'  for details about URL formats --
#'  including local file addresses -- that might work here.
#'  
#'  \code{\link{import_biom}}
#'  
#'  \code{\link{import_qiime}}
#'  
#' @export
#' @examples
#' # This should return TRUE on your system if you have internet turned on
#' # and a standard R installation. Indicates whether this is likely to
#' # work on your system for a URL or local file, respectively.
#' capabilities("http/ftp"); capabilities("fifo")
#' # A working example with a local example file included in phyloseq
#' zipfile = "study_816_split_library_seqs_and_mapping.zip"
#' zipfile = system.file("extdata", zipfile, package="phyloseq")
#' tarfile = "study_816_split_library_seqs_and_mapping.tar.gz"
#' tarfile = system.file("extdata", tarfile, package="phyloseq")
#' tarps = microbio_me_qiime(tarfile)
#' zipps = microbio_me_qiime(zipfile) 
#' identical(tarps, zipps)
#' tarps; zipps
#' plot_heatmap(tarps)
#' # A real example 
#' # # Smokers dataset
#' # smokezip = "ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_524_split_library_seqs_and_mapping.zip"
#' # smokers1 = microbio_me_qiime(smokezip)
#' # # Alternatively, just use the study number
#' # smokers2 = microbio_me_qiime(524)
#' # identical(smokers1, smokers2)
microbio_me_qiime = function(zipftp, ext=".zip", parsef=parse_taxonomy_greengenes, ...){
	# Define naming convention
	front = "ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_"	
	if( inherits(zipftp, "numeric") ){
		# If study number instead of string,
		# create the ftp URL using ext and convention
		back  = paste0("_split_library_seqs_and_mapping", ext)
		zipftp = paste0(front, zipftp, back)
	} else {
		# Determine file extension from the file path itself
		ext = substring(zipftp, regexpr("\\.([[:alnum:]]+)$", zipftp)[1])
		back  = paste0("_split_library_seqs_and_mapping", ext)
	}
	# Check if zipftp is clearly an externally located file, ftp, http, etc.
	externprefixes = c("http://", "https://", "ftp://")
	prefix = regexpr("^([[:alnum:]]+)\\://", zipftp)
	if( substr(zipftp, 1, attr(prefix, "match.length")[1]) %in% externprefixes ){
		# If external, then create temporary file and download
		zipfile = tempfile()
		download.file(zipftp, zipfile, ...)
	} else {
		# Else it is a local zipfile
		zipfile = zipftp
	}
	# Use the apparent file naming convention for microbio.me/qiime
	# as the de facto guide for this API. In particular,
  # the expectation o fthe study name (already used above)
	studyname = gsub("\\_split\\_.+$", "", basename(zipftp))
  # The output of tempdir() is always the same in the same R session
  # To avoid conflict with multiple microbio.me/qiime unpacks 
  # in the same session, pre-pend the study name and datestamp
  unpackdir = paste0(studyname, "_", gsub("[[:blank:][:punct:]]", "", date()))
  # Add the temp path
	unpackdir = file.path(tempdir(), unpackdir)
  # Create the unpack directory if needed (most likely).
  if( !file.exists(unpackdir) ){dir.create(unpackdir)}
	# Unpack to the temporary directory using unzip or untar
	if( ext == ".zip" ){	
		unzip(zipfile, exdir=unpackdir, overwrite=TRUE)
	} else if( ext %in% c("tar.gz", ".tgz", ".gz", ".gzip", ".bzip2", ".xz") ){
		# untar the tarfile to the new temp dir
		untar(zipfile, exdir=unpackdir)
	} else {
		# The compression format was not recognized. Provide informative error msg.
    msg = paste("Could not determine the compression type.",
                "Expected extensions are (mostly):",
                ".zip, .tgz, .tar.gz", sep="\n")
		stop(msg)
	}  
  # Define a list of imported objects that might grow 
  # if the right file types are present and imported correctly.
	imported_objects = vector("list")
	# Search recursively in the unpacked directory for the .biom file
  # and parse if it is.
	# There should be only one. Throw warning if more than one, take the first.
	biomfile = list.files(unpackdir, "\\.biom", full.names=TRUE, recursive=TRUE)
  if( length(biomfile) > 1 ){
    warning("more than one .biom file found in compressed archive. Importing first only.")
    biomfile = biomfile[1]
  } else if( length(biomfile) == 1 ){
		cat("Found biom-format file, now parsing it... \n")
		biom = import_biom(biomfile, parseFunction=parsef)
		cat("Done parsing biom... \n")
		imported_objects = c(imported_objects, list(biom))
	}
	# Check if sample_data (qiime mapping) file present, and parse if it is.
	sdfile = list.files(unpackdir, "\\_mapping\\_file\\.txt", full.names=TRUE, recursive=TRUE)
	if( length(sdfile) > 1 ){
	  warning("more than one mapping file found in compressed archive. Importing first only.")
	  sdfile = sdfile[1]
	} else if( length(sdfile)==1 ){
    cat("Importing Sample Metdadata from mapping file...", fill=TRUE)
		sample_metadata = import_qiime_sample_data(sdfile)
		imported_objects = c(imported_objects, list(sample_metadata))
	}
	# Check success, notify user, and return.
	if( length(imported_objects) > 1 ){
		# If there are more than one imported objects, merge them and return
		cat("Merging the imported objects... \n")
		physeq = do.call("merge_phyloseq", imported_objects)
		if( inherits(physeq, "phyloseq") ){
			cat("Successfully merged, phyloseq-class created. \n Returning... \n")	
		}
		return(physeq)
	} else if( length(imported_objects) == 1 ){
		cat("Note: only on object in the zip file was imported. \n")
		cat("It was ", class(imported_objects[[1]]), " class. \n")
		return(imported_objects[[1]])
	} else {
		cat("PLEASE NOTE: No objects were imported. \n", 
				"You chould check the zip file, \n",
				"as well as the naming conventions in the zipfile \n",
				"to make sure that they match microbio.me/qiime. \n",
				"Instead returning NULL... \n")
		return(NULL)
	}
}
################################################################################
#' Import usearch table format (\code{.uc}) to OTU table
#' 
#' UPARSE is an algorithm for OTU-clustering implemented within usearch.
#' At last check, the UPARSE algortihm was accessed via the 
#' \code{-cluster_otu} option flag.
#' For details about installing and running usearch, please refer to the
#' \href{http://drive5.com/usearch/}{usearch website}.
#' This importer is intended to read a particular table format output
#' that is generated by usearch with the \code{-uc} option flag,
#' a file format that is often given the \code{.uc} extension
#' in usearch documentation.
#' Because usearch is an external (non-R) application, there is no direct
#' way to continuously check that these suggested arguments and file formats will 
#' remain in their current state. 
#' If there is a problem, please verify your version of usearch,
#' create a small reproducible example of the problem,
#' and post it as an issue on the phyloseq issues tracker.
#' The version of usearch expected by this import function is \code{7.0.109}.
#' Hopefully later versions of usearch maintain this function and format,
#' but the phyloseq team has no way to guarantee this,
#' and so any feedback about this will help maintain future functionality.
#' For instance, it is currently
#' assumed that the 9th and 10th columns of the \code{.uc} table
#' hold the read-label and OTU ID, respectively;
#' and it is also assumed that the delimiter between sample-name and read
#' in the read-name entries is a single \code{"_"}.
#' 
#' @param ucfile (Required). A file location character string 
#'  or \code{\link{connection}}
#'  corresponding to the file that contains the usearch output table.
#'  This is passed directly to \code{\link{read.table}}.
#'  Please see its \code{file} argument documentation for further
#'  links and details.
#' 
#' @param colRead (Optional). Numeric. The column index in the uc-table 
#'  file that holds the read IDs.
#'  The default column index is \code{9}.
#' 
#' @param colOTU (Optional). Numeric. The column index in the uc-table 
#'  file that holds OTU IDs.
#'  The default column index is \code{10}.
#'   
#' @param readDelimiter (Optional). An R \code{\link{regex}} as a character string.
#'  This should be the delimiter that separates the sample ID
#'  from the original ID in the demultiplexed read ID of your sequence file.
#'  The default is plain underscore, which in this \code{\link{regex}} context
#'  is \code{"_"}.
#'  
#' @param verbose (Optional). A \code{\link{logical}}.
#'  Default is \code{TRUE}. 
#'  Should progresss messages
#'  be \code{\link{cat}}ted to standard out?
#' 
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @export
#' @seealso \code{\link{import}}
#' 
#' \code{\link{import_biom}}
#' 
#' \code{\link{import_qiime}}
#' 
#' @examples
#' usearchfile <- system.file("extdata", "usearch.uc.gz", package="phyloseq")
#' import_usearch_uc(usearchfile)
import_usearch_uc <- function(ucfile, colRead=9, colOTU=10,
                              readDelimiter="_", verbose=TRUE){
  if(verbose){cat("Reading `ucfile` into memory and parsing into table \n")}
  x = fread(paste0(readLines(ucfile), collapse="\n"), 
            sep="\t", header=FALSE, colClasses="character", 
            na.strings=c("*", '*', "NA","N/A",""))
  colNames = colnames(x)
  colNames[c(colRead, colOTU)] <- c("read", "OTU")
  setnames(x, colNames)
  NrawEntries = nrow(x)
  if(verbose){
    cat("Initially read", NrawEntries, "entries. \n")
    cat("... Now removing unassigned OTUs (* or NA)... \n")
  }
  #x = subset(x, !is.na(OTU))
  x = x[!is.na(x$OTU), ]
  if(verbose){
    cat("Removed", NrawEntries - nrow(x), "entries that had no OTU assignment. \n")
    cat("A total of", nrow(x), "will be assigned to the OTU table.\n")
  }
  # Process sequence label to be sample label only
  x$read <- gsub(paste0(readDelimiter, ".+$"), "", x$read) 
  # Convert long (melted) table into a sample-by-OTU OTU table, and return
  return(otu_table(as(table(x$read, x$OTU), "matrix"), taxa_are_rows=FALSE))
}
################################################################################
################################################################################
################################################################################
################################################################################