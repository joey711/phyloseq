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
#' @param ... (Required). Additional arguments providing file paths, and possible
#'  other paramaters to the desired tool-specific import function.
#'
#' @return In most cases a \code{\link{phyloseq-class}} will be returned, though
#'  the included component data will vary by pipeline/tool, and also
#'  by the types of data files provided.
#'  The expected behavior is to return the most-comprehensive object possible,
#'  given the provided arguments and pipeline/tool.
#'
#' @seealso 
#' For mothur, see:
#' \code{\link{import_mothur}}
#' 
#' Separate tools for mothur are also:
#' \code{\link{show_mothur_list_cutoffs}}
#' \code{\link{import_mothur_dist}}
#' \code{\link{export_mothur_dist}}
#' 
#' For PyroTagger, see:
#' \code{\link{import_pyrotagger_tab}}
#'
#' For QIIME, see:
#' \code{\link{import_qiime}}
#'
#' For BIOM format, see:
#' \code{\link{import_biom}}
#' 
#' For RDP pipeline, see:
#' \code{\link{import_RDP_cluster}}
#'
#' \code{\link{import_RDP_otu}}
#' 
#' @references 
#' mothur: \url{http://www.mothur.org/wiki/Main_Page}
#'
#' PyroTagger: \url{http://pyrotagger.jgi-psf.org/}
#'
#' QIIME: \url{http://qiime.org/}
#'
#' BIOM: \url{http://www.biom-format.org/}
#'
#' RDP pipeline: \url{http://pyro.cme.msu.edu/index.jsp}
#' 
#' @export
#' @examples
#'  ## import("QIIME", otufilename=myOtuTaxFilePath, mapfilename=myMapFilePath)
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
#' Import function to read files created by the QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' including especially an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information. The map-file is also an important input
#' to QIIME that stores sample covariates, converted naturally to the 
#' \code{\link{sample_data-class}} component data type in the phyloseq-package. 
#' QIIME may also produce a
#' phylogenetic tree with a tip for each OTU, which can also be imported by this
#' function.
#' 
#' See \url{"http://www.qiime.org/"} for details on using QIIME. While there are
#' many complex dependencies, QIIME can be downloaded as a pre-installed 
#' linux virtual machine that runs ``off the shelf''. 
#'
#' The different files useful for import to \emph{phyloseq} are not collocated in
#' a typical run of the QIIME pipeline. See the main \emph{phyloseq} vignette for an 
#' example of where ot find the relevant files in the output directory. 
#'
#' @usage import_qiime(otufilename=NULL, mapfilename=NULL,
#'  treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
#'  parseFunction=parse_taxonomy_qiime, showProgress=TRUE, chunk.size=1000L, ...)
#'
#' @param otufilename (Optional). A character string indicating the file location of the OTU file.
#' The combined OTU abundance and taxonomic identification file,
#' tab-delimited, as produced by QIIME under default output settings.
#'  Default value is \code{NULL}. 
#' 
#' @param mapfilename (Optional). The QIIME map file required for processing pyrosequencing tags
#' in QIIME as well as some of the post-clustering analysis. This is a required
#' input file for running QIIME. Its strict formatting specification should be
#' followed for correct parsing by this function.
#'  Default value is \code{NULL}. 
#'
#' @param treefilename (Optional).
#'  A file representing a phylogenetic tree
#'  or a \code{\link{phylo}} object.
#'  Files can be NEXUS or Newick format.
#'  See \code{\link{read_tree}} for more details. 
#'  If provided, the tree should have the same OTUs/tip-labels
#'  as the OTUs in the other files. 
#'  Any taxa or samples missing in one of the files is removed from all. 
#'  For the QIIME pipeline
#'  this tree is typically a tree of the representative 16S rRNA sequences from each OTU
#'  cluster, with the number of leaves/tips equal to the number of taxa/species/OTUs.
#'  Default value is \code{NULL}. 
#'  Note that this argument can be a tree object (\code{\link[ape]{phylo}}-class)
#'  for cases where the tree has been --- or needs to be --- imported separately.
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
#' @param showProgress (Optional). A logical. 
#'  Indicates whether import progress/status should be printed to the terminal.
#'  Default value is \code{TRUE}, meaning the progress will be shown. 
#'
#' @param chunk.size (Optional). Positive integer. Default is \code{1000L}.
#'  This is the number of lines to be read and processed at a time,
#'  passed directly to the \code{\link{import_qiime_otu_tax}} function.
#'  A lower value helps control memory-fault, but slows down the import.
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
#' @import Biostrings
#' @export
#' @examples
#'  otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
#'  mapfile <- system.file("extdata", "master_map.txt", package="phyloseq")
#'  trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
#'  import_qiime(otufile, mapfile, trefile, showProgress=FALSE)
import_qiime <- function(otufilename=NULL, mapfilename=NULL,
	treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
	parseFunction=parse_taxonomy_qiime, showProgress=TRUE, chunk.size=1000L, ...){

	# initialize the argument-list for phyloseq. Start empty.
	argumentlist <- list()

	if( !is.null(mapfilename) ){	
		if( showProgress==TRUE ){
			cat("Processing map file...", fill=TRUE)
		}
		QiimeMap     <- import_qiime_sample_data(mapfilename)
		argumentlist <- c(argumentlist, list(QiimeMap))
	}

	if( !is.null(otufilename) ){
		if( showProgress==TRUE ){
			cat("Processing otu/tax file...", fill=TRUE)
		}		
		otutax <- import_qiime_otu_tax(otufilename, parseFunction, FALSE, chunk.size, showProgress)
		otutab <- otu_table(otutax$otutab, TRUE)
		tax_table <- tax_table(otutax$tax_table)
		argumentlist <- c(argumentlist, list(otutab), list(tax_table) )
	}

	if( !is.null(treefilename) ){
		if( showProgress==TRUE ){
			cat("Processing phylogenetic tree...", fill=TRUE)
		}
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
	
	if( !is.null(refseqfilename) ){
		if( showProgress==TRUE ){
			cat("Processing Reference Sequences...", fill=TRUE)
		}
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
#' read_tree(treefile, errorIfNULL=FALSE, ...)
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
#'  be hard to track in treefileyour own code if you're not aware of this
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
#'	\code{\link{phylo}}, \code{\link[ape]{read.tree}}, \code{\link[ape]{read.nexus}}
#'
#' @import ape
#' @export
#' @examples
#' read_tree(system.file("extdata", "esophagus.tree.gz", package="phyloseq"))
#' read_tree(system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq"))
read_tree <- function(treefile, errorIfNULL=FALSE, ...){
	# "phylo" object provided directly
	if( class(treefile)[1] %in% c("phylo") ){ 
		tree <- treefile
	# file path to tree file provided (NEXUS)
	} else {
		# # Try Nexus first, protected, then newick if it fails
		tree <- NULL
		try(tree <- read.nexus(treefile, ...), TRUE)
		# Try Newick if nexus didn't work.
		if(is.null(tree)) try(tree <- read.tree(treefile, ...), TRUE)
	}
	# If neither tree-import worked (still NULL), report warning
	if( errorIfNULL & is.null(tree) ){
		stop("tree file could not be read.\nPlease retry with valid tree.")
	}
	return(tree)
}
################################################################################
#' Import a QIIME-formatted otu-tax file into a list of two matrices.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' including especially an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information.   
#' 
#' See \url{"http://www.qiime.org/"} for details on using QIIME. While there are
#' many complex dependencies, QIIME can be downloaded as a pre-installed 
#' linux virtual machine that runs ``off the shelf''. 
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
#' @usage import_qiime_otu_tax(file, parseFunction=parse_taxonomy_qiime, parallel=FALSE, chunk.size=1000, verbose=TRUE)
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
#' @param parallel (Optional). Logical. Should the parsing be performed in 
#'  parallel?. Default is \code{FALSE}. Only a few steps are actually 
#'  parallelized, and for most datasets it will actually be faster and 
#'  more efficient to keep this set to \code{FALSE}.
#'  Also, to get any benefit at all, you will need to register a 
#'  parallel ``backend'' through one of the backend packages supported
#'  by the \code{\link{foreach-package}}.
#'
#' @param chunk.size (Optional). Positive integer. Default is \code{1000L}.
#'  This is the number of lines to be read and processed at a time. In many
#'  cases this will have a much larger effect on speed and efficiency than
#'  the setting for \code{parallel}. If you set this value too large, you
#'  risk a memory fault. If you set it too low, you risk a very slow
#'  parsing process.
#'
#' @param verbose (Optional). Logical. Default is \code{TRUE}. 
#'  Should status updates of the parsing process be printed to screen?
#'
#' @return A list of two matrices. \code{$otutab} contains the OTU Table
#'  as a numeric matrix, while \code{$tax_table} contains a character matrix
#'  of the taxonomy assignments.
#'
#' @import foreach
#' @importFrom plyr llply
#'
#' @seealso \code{\link{import}}, \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
#'  \code{\link{import_qiime}}
#'  \code{\link{import_qiime_otu_tax}}
#'  \code{\link{import_env_file}}
#' @export
#' @examples
#'  otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
#'  import_qiime_otu_tax(otufile)
import_qiime_otu_tax <- function(file, parseFunction=parse_taxonomy_qiime, parallel=FALSE, chunk.size=1000L, verbose=TRUE){

	### Some parallel-foreach housekeeping.    
    # If user specifies not-parallel run (the default), register the sequential "back-end"
    if( !parallel ){ registerDoSEQ() }

	# Check for commented lines, starting with line 1.
	# The deepest commented line (largest n) is assumed to have header information.
	n <- 1L
	header      <- readLines(file, n)
	stillHeader <- identical(substr(header[1], 1, 1), "#")
	while( stillHeader ){
		# Check again if this is a commented line
		n <- n + 1
		header      <- strsplit(readLines(file, n)[[n]], "\t", TRUE)[[1]]
		stillHeader <- identical(substr(header[1], 1, 1), "#")
		# Read the header portion of file up to the line
		# that contains the "Consensus Lineage" label
		stillHeader <- stillHeader & !any(c("Consensus Lineage", "Consensus", "Lineage") %in% header)
	}

	# Remove the first and last values (strings "#OTU ID" and "Consensus Lineage", respectively)
	# to get the sample names from the header
	sampleNames <- header[-length(header)][-1]

	if(verbose) cat("\nReading and parsing file in chunks ... Could take some time. Please be patient...", fill=TRUE)
	if(verbose) cat("\nBuilding OTU Table in chunks. Each chunk is one dot.", fill=TRUE)
	
	# Initialize taxstring, you will just build it up and parse it at once after the loop.
	taxstring <- character(0)
	
	# Initialize the empty matrix to build-up as you go.
	otutab <- matrix(NA_integer_, nrow=0, ncol=length(sampleNames))
	
	# Loop while chunk is as big as chunk.size
	chunking_otu <- TRUE
	while( chunking_otu ){
		# optionally print chunk indicator
		if(verbose) cat(".")
			
		# Define the chunk
		taxa.scan <- scan(file, as.list(header), skip=n, quiet=TRUE, sep="\t", multi.line=FALSE, nmax=chunk.size)
		
		# Store the taxa-names and taxonomy-string separately.
		taxaNames <- taxa.scan[[1]]
		taxstring <- c(taxstring, taxa.scan[[length(taxa.scan)]])

		# Trim taxa.scan to just abundance values, and name by sample
		taxa.scan <- taxa.scan[-length(header)][-1]
		names(taxa.scan) <- sampleNames	
		
		# Parse otu_table. 
		otutab.chunk <- foreach( i=sampleNames, .combine=cbind) %dopar%  {as.numeric(taxa.scan[[i]])}
		colnames(otutab.chunk) <- sampleNames 
		rownames(otutab.chunk) <- taxaNames	
		
		# Add to otutab
		otutab <- rbind(otutab, otutab.chunk)
		
		# Control loop. Either stop, or calculate new n.
		if( length(taxa.scan[[1]]) < chunk.size ){ 
			chunking_otu <- FALSE
		} else {
			# Calculate the new start position, n
			n <- n + chunk.size
		}
	}

	# Remove taxa.scan to clear memory usage. Call garbage collection right away.
	rm(taxa.scan)
	garbage.collection <- gc(FALSE)
	
	if(verbose) cat("Building Taxonomy Table...", fill=TRUE)
	
	# # # # # Explicit foreach loops lost the competition (by a lot) to llply
	# Split into "jagged" list (vectors of different lengths)
	taxlist = llply(taxstring, parseFunction, .parallel=parallel)
	# Add OTU names to list element names
	names(taxlist) <- rownames(otutab)
	# Build the tax table from the jagged list.	
	tax_table <- build_tax_table(taxlist)

	# Call garbage collection one more time. Lots of unneeded stuff.
	garbage.collection <- gc(FALSE)
	return(list(otutab=otutab, tax_table=tax_table))

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
#' @seealso \code{\link{import}}, \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
#'  \code{\link{import_qiime}}
#'  \code{\link{import_qiime_otu_tax}}
#'  \code{\link{import_env_file}}
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
######################################################################################
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
#' @seealso \code{\link{import}}, \code{\link{tip_glom}}
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
#' @seealso \code{\link{import_env_file}}, \code{\link{tip_glom}}, \code{\link{otu_table}}
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
################################################################################
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
################################################################################
################################################################################
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
#' @seealso \code{\link{import_mothur}}
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
#' @keywords internal
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
	# cutoff is often numeric, but the index itself must be a character.
	if( class(cutoff)=="numeric" ){ cutoff <- as(cutoff, "character") }
	
	# The first two elements are the cutoff and the number of OTUs. metadata. rm.
	OTUs <- tabsplit[[cutoff]][-(1:2)]
	
	# split each element on commas
	OTUs <- strsplit(OTUs, ",", fixed=TRUE)
	
	# Name each OTU (currently as the first seq name in each cluster), and return the list
	# # # names(OTUs) <- paste("OTUID_", 1:length(OTUs), sep="")
	names(OTUs) <- sapply(OTUs, function(i){return(i[1])})
	
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
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
#'
import_mothur_groups <- function(mothur_group_file){
	group_table <- read.table(mothur_group_file, sep="\t")
	rownames(group_table) <- as(group_table[, 1], "character")
	return(group_table)
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
#'  by the \code{mothur_list_file} argument). The default
#'  is to take the largest value among the cutoff values contained in the list
#'  file. If only one cutoff is included in the file, it is taken and this
#'  argument does not need to be specified. Note that the \code{cluster()}
#'  function within the \emph{mothur} package will often produce a list file
#'  with multiple cutoff values, even if a specific cutoff is specified. It is
#'  suggested that you check which cutoff values are available in a given list
#'  file using the \code{\link{show_mothur_list_cutoffs}} function.
#'
#' @return An \code{\link{otu_table}} object.
#'
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
#'  
import_mothur_otu_table <- function(mothur_list_file, mothur_group_file, cutoff=NULL){
	
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
	
	# Finally, return the otu_table as a phyloseq otu_table object.
	return(otu_table(mothur_otu_table, taxa_are_rows=TRUE))
}
################################################################################
#' Import and prune mothur-produced tree.
#'
#' The \code{\link[ape]{read.tree}} function is sufficient for importing a 
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
#' The \code{import_mothur} function
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
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
#' 
import_mothur_tree <- function(mothur_tree_file, mothur_list_file, cutoff=NULL){
	
	otulist <- import_mothur_otulist(mothur_list_file, cutoff)

	# Read the original all-sequences tree from mothur
	tree <- read_tree(mothur_tree_file)
	
	# Loop to merge the sequences in the tree by OTU
	# cycle through each otu, and sum the number of seqs observed for each sample
	for( i in names(otulist) ){
		# i <- names(otulist)[1]
		# First merge the reads that are in the same OTU ("eqtaxa" argument)
		tree <- merge_taxa(tree, otulist[[i]])
		# Rename the tip that was kept to the otuID.
		# By default, the first element of eqtaxa is used as archetype.
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
#'  species/taxa abundance table (\code{otu_table}). See
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
#'  \code{otu_table}-class abundance table. If only a list and group file are 
#'  provided, then an \code{otu_table} object is returned. Similarly, if only a list
#'  and tree file are provided, then only a tree is returned (\code{"phylo"} class).
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
#' # show_mothur_list_cutoffs(mothur_list_file)
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
import_mothur <- function(mothur_list_file,  mothur_group_file=NULL,
		mothur_tree_file=NULL, cutoff=NULL){

	if( missing(mothur_list_file) ){cat("you must provide the mothur_list_file argument\n")}

	# Create otu_table object, OTU, only if group file provided. 
	if( !is.null(mothur_group_file) ){
		OTU  <- import_mothur_otu_table(mothur_list_file, mothur_group_file, cutoff)		
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
#' @examples #
#' ### data(GlobalPatterns) 
#' ### myDistObject <- as.dist(cophenetic(tre(GlobalPatterns)))
#' ### export_mothur_dist(myDistObject, "myfilepathname.dist")
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
#' @import ape
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
#' @param treefilename (Optional).
#'  A file representing a phylogenetic tree
#'  or a \code{\link{phylo}} object.
#'  Files can be NEXUS or Newick format.
#'  See \code{\link{read_tree}} for more details. 
#'  If provided, the tree should have the same OTUs/tip-labels
#'  as the OTUs in the other files. 
#'  Any taxa or samples missing in one of the files is removed from all. 
#'  For the QIIME pipeline
#'  this tree is typically a tree of the representative 16S rRNA sequences from each OTU
#'  cluster, with the number of leaves/tips equal to the number of taxa/species/OTUs.
#'  Default value is \code{NULL}. 
#'  Note that this argument can be a tree object (\code{\link[ape]{phylo}}-class)
#'  for cases where the tree has been --- or needs to be --- imported separately.
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
#'  As the BIOM format evolves, version-specific importers will be available
#'  by adjusting the version value. Default is \code{0.9}. Not implemented.
#'  Has no effect (yet).
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
#' \code{\link[Biostrings]{XStringSet-io}}
#'
#' @references \href{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{biom-format}
#'
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr ldply
#' @importFrom plyr laply
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
	x = fromJSON(BIOMfilename)
	
	########################################
	# OTU table:
	########################################
	# Check if sparse. Must parse differently than dense
	if( x$matrix_type == "sparse" ){
		otumat = matrix(0, nrow=x$shape[1], ncol=x$shape[2])
		dummy  = sapply(x$data, function(i){otumat[(i[1]+1), (i[2]+1)] <<- i[3]})
	}
	# parse the dense matrix instead.
	if( x$matrix_type == "dense" ){
		# each row will be complete data values, should use laply
		otumat = laply(x$data, function(i){i}, .parallel=parallel)
	}
	
	# Get row (OTU) and col (sample) names
	rownames(otumat) = sapply(x$rows, function(i){i$id})
	colnames(otumat) = sapply(x$columns, function(i){i$id})
	
	otutab = otu_table(otumat, TRUE)
	
	argumentlist <- c(argumentlist, list(otutab))
	
	########################################
	# Taxonomy Table
	########################################
	# Need to check if taxonomy information is empty (minimal BIOM file)
	if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
		tax_table <- NULL
	} else {
		# parse once each character vector, save as a list
		taxlist = lapply(x$rows, function(i){
			parseFunction(i$metadata$taxonomy)
		})
		names(taxlist) = sapply(x$rows, function(i){i$id})
		tax_table = build_tax_table(taxlist)
	}
	argumentlist <- c(argumentlist, list(tax_table))
	
	########################################
	# Sample Data ("columns" in QIIME/BIOM)
	########################################
	# If there is no metadata (all NULL), then set sam_data <- NULL
	if(  all( sapply(sapply(x$columns, function(i){i$metadata}), is.null) )  ){
		sam_data <- NULL
	} else {
		sam_data <- ldply(x$columns, function(i){
			if( class(i$metadata) == "list"){
				return(i$metadata[[1]])
			} else {
				return(i$metadata)				
			}
		}, .parallel=parallel)
		rownames(sam_data) <- sapply(x$columns, function(i){i$id})
		sam_data <- sample_data(sam_data)
	}
	argumentlist <- c(argumentlist, list(sam_data))

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