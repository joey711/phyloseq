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
#' For RDP pipeline, see:
#' \code{\link{import_RDP_cluster}}
#' 
#' @references 
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
######################################################################################
#' Import function to read files created by the QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' including especially an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information. The map-file is also an important input
#' to QIIME that stores sample covariates, converted naturally to the 
#' \code{\link{sampleData-class}} component data type in the phyloseq-package. 
#' QIIME may also produce a
#' phylogenetic tree with a tip for each OTU, which can also be imported by this
#' function.
#' 
#' See \code{"http://www.qiime.org/"} for details on using QIIME. While there are
#' many complex dependencies, QIIME can be downloaded as a pre-installed 
#' linux virtual machine that runs ``off the shelf''. 
#'
#' The different files useful for import to \emph{phyloseq} are not collocated in
#' a typical run of the QIIME pipeline. See the main \emph{phyloseq} vignette for an 
#' example of where ot find the relevant files in the output directory. 
#'
#' @usage import_qiime(otufilename=NULL, mapfilename=NULL,
#'	treefilename=NULL, biotaxonomy=NULL, ...)
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
#' @param treefilename (Optional). A phylogenetic tree in NEXUS format. For the QIIME pipeline
#'  this is typically a tree of the representative 16S rRNA sequences from each OTU
#'  cluster, with the number of leaves/tips equal to the number of taxa/species/OTUs.
#'  Default value is \code{NULL}. ALTERNATIVELY, this argument can be a tree object
#'  (\code{\link[ape]{phylo}}-class), in case the tree has already been
#'  imported, or is in a different format than NEXUS.
#'
#' @param biotaxonomy (Optional). A character vector indicating the name of each taxonomic level
#'  in the taxonomy-portion of the otu-file, which may not specify these levels 
#'  explicitly.
#'  Default value is \code{NULL}. 
#'
#' @param ... Additional arguments passed to \code{\link[ape]{read.nexus}}, as necessary.
#'  Make sure that your phylogenetic tree file is readable by \code{\link[ape]{read.nexus}}
#'  prior to calling this function.
#'
#' @return A \code{\link{phyloseq-class}} object.
#'
#' @seealso \code{\link{phyloseq}}, \code{\link{merge_phyloseq}},
#'	 \code{\link[ape]{read.tree}}, \code{\link[ape]{read.nexus}}
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
#' @export
#' @examples
#'  # import_qiime(myOtuTaxFilePath, myMapFilePath)
import_qiime <- function(otufilename=NULL, mapfilename=NULL,
	treefilename=NULL, biotaxonomy=NULL, ...){

	# initialize the argument-list for phyloseq. Start empty.
	argumentlist <- list()
		
	if( is.null(biotaxonomy) ){
	 	biotaxonomy=c("Root", "Domain", "Phylum", "Class", "Order",
		 	"Family", "Genus", "Species", "Strain")
	 }

	if( !is.null(mapfilename) ){	
		# Process mapfile. Name rows as samples.
		QiimeMap     <- import_qiime_sampleData(mapfilename)
		argumentlist <- c(argumentlist, list(QiimeMap))
	}

	if( !is.null(otufilename) ){
		otutax       <- import_qiime_otu_tax(otufilename, biotaxonomy)	
		argumentlist <- c(argumentlist, list(otuTable(otutax)), list(taxTab(otutax)) )
	}

	if( !is.null(treefilename) ){
		if( class(treefilename) %in% c("phylo") ){ 
			tree <- treefilename
		} else {
			#tree <- ape::read.tree(treefilename, ...)
			tree <- ape::read.nexus(treefilename, ...)
		}
		argumentlist <- c(argumentlist, list(tree) )
	}

	do.call("phyloseq", argumentlist)
}
######################################################################################
#' Import just OTU/Taxonomy file from QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' including especially an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information.   
#' 
#' Add reference to the QIIME pipeline.
#'
#' @param otufilename A character string indicating the file location of the OTU file.
#' The combined OTU abundance and taxonomic identification file,
#' tab-delimited, as produced by QIIME under default output settings.
#'  Default value is \code{NULL}. 
#' 
#' @param biotaxonomy A character vector indicating the name of each taxonomic level
#'  in the taxonomy-portion of the otu-file, which may not specify these levels 
#'  explicitly.
#'  Default value is \code{NULL}. 
#'
#' @return An \code{otuTax} object.
#'
#' @seealso \code{\link{merge_phyloseq}}, \code{\link{phyloseq}}, 
#'   \code{\link{import_qiime_sampleData}}
#'
#' @keywords internal
import_qiime_otu_tax <- function(otufilename, biotaxonomy=NULL){
	if( is.null(biotaxonomy) ){
	 	biotaxonomy=c("Root", "Domain", "Phylum", "Class", "Order",
		 	"Family", "Genus", "Species", "Strain")
	 }
	##########################################
	# Process otu table. "otuID" convention
	# specific to Qiime. Might need to be abstracted.
	##########################################
	# first read the table. Skip line 1, avoid comment character "#"
	otutab <- read.table(file=otufilename, header=TRUE,
		sep="\t", comment.char="", skip=1)
	# name the rows by otuID
	rownames(otutab) <- paste("otuID", as.character(otutab[,1]), sep="_")
	# remove the otuID column
	otutab <- otutab[, 2:ncol(otutab)]
	##########################################
	# Process/separate the lineage information
	##########################################
	splitaxa <- strsplit(as.character(otutab[,"Consensus.Lineage"]),";")
	taxtab   <- matrix(NA, nrow(otutab), length(biotaxonomy))
	colnames(taxtab) <- biotaxonomy
	# If present, place the taxonomy labels in matrix
	# starting on the left column (root)
	for( i in 1:nrow(otutab) ){
		taxtab[i, 1:length(splitaxa[[i]])] <- splitaxa[[i]]
	}
	rownames(taxtab) <- rownames(otutab)
	taxtab <- taxTab( as.matrix(taxtab) )
	
	# Remove taxonomy column from otutab
	otutab <- otutab[, which(colnames(otutab)!="Consensus.Lineage")]
	# convert to matrix of integers, and then otuTable object
	otutab <- otuTable( as(otutab, "matrix"), speciesAreRows=TRUE )

	return( phyloseq(otutab, taxtab) )
}
######################################################################################
#' Import just \code{sampleData} file from QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' This includes the map-file, which is an important \emph{input}
#' to QIIME that can also indicate sample covariates. It is converted naturally to the 
#' sampleData component data type in phyloseq-package, based on the R data.frame.   
#' 
#' See \code{\link{import_qiime}} for more information about QIIME. It is also the
#' suggested function for importing QIIME-produced data files. 
#'
#' @usage import_qiime_sampleData(mapfilename)
#'
#' @param mapfilename (Required). A character string. The name of the QIIME map
#'  file required for processing pyrosequencing tags
#'  in QIIME as well as some of the post-clustering analysis. This is a required
#'  input file for running QIIME. Its strict formatting specification is expected by
#'  this function, do not attempt to modify it manually once it has worked properly
#'  in QIIME. 
#'
#' @return A \code{sampleData} object.
#'
#' @seealso \code{\link{import_qiime}}, \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
#'	 \code{\link{import_qiime_otu_tax}}
#'
#' @keywords internal
import_qiime_sampleData <- function(mapfilename){
	# Process mapfile. Name rows as samples.
	QiimeMap <- read.table(file=mapfilename, header=TRUE,
		sep="\t", comment.char="")
	rownames(QiimeMap) <- as.character(QiimeMap[,1])
	return( sampleData(QiimeMap) )
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
#'  the output otuTable. 
#'
#' @param sep A character string indicating the delimiter used in the file.
#'  The default is \code{"\t"}.
#'
#' @param ... Additional parameters passed on to \code{\link{read.table}}.
#'
#' @return An \code{\link{otuTable-class}}, or \code{\link{phyloseq-class}} if 
#'  a \code{\link{phylo-class}} argument is provided to \code{tree}.
#'
#' @references \url{http://bmf2.colorado.edu/unifrac/}
#' 
#' @seealso \code{\link{import}}, \code{\link{tipglom}}
#' @export
#' @examples 
#' # import_env_file(myEnvFile, myTree)
import_env_file <- function(envfilename, tree=NULL, sep="\t", ...){
	tipSampleTable <- read.table(envfilename, sep=sep, ...)
	# Convert to otuTable-class table (trivial table)
	physeq <- envHash2otuTable(tipSampleTable)
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
#' table as arguments to \code{\link{tipglom}} to create an object with a
#' non-trivial \code{otuTable}.  
#'
#' @usage envHash2otuTable(tipSampleTable)
#'
#' @param tipSampleTable (Required). A two-column character table (matrix or data.frame), 
#' where each row specifies the sequence name and source sample, consistent with the 
#' env-file for the UniFrac server (\url{http://bmf2.colorado.edu/unifrac/}). 
#'
#' @return \code{\link{otuTable}}. A trivial OTU table where each sequence 
#'  is treated as a separate OTU. 
#' 
#' @references \url{http://bmf2.colorado.edu/unifrac/}
#' 
#' @seealso \code{\link{import_env_file}}, \code{\link{tipglom}}, \code{\link{otuTable}}
#'
#' @keywords internal
#' @examples #
#' ## fakeSeqNameVec <- paste("seq_", 1:8, sep="")
#' ## fakeSamNameVec <- c(rep("A", 4), rep("B", 4))
#' ## fakeSeqAbunVec <- sample(1:50, 8, TRUE)
#' ## test    <- cbind(fakeSeqNameVec, fakeSamNameVec, fakeSeqAbunVec)
#' ## testotu <- envHash2otuTable( test )
#' ## test    <- cbind(fakeSeqNameVec, fakeSamNameVec)
#' ## testotu <- envHash2otuTable( test )
envHash2otuTable <- function(tipSampleTable){
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
	return( otuTable(trivialOTU, speciesAreRows=TRUE) )
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#' Import RDP cluster file and return otuTable (abundance table).
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
#' coerced to an \code{\link{otuTable}} object, and returned.
#'
#' @usage import_RDP_cluster(RDP_cluster_file)
#'
#' @param RDP_cluster_file A character string. The name of the \code{".clust"} 
#'  file produced by the 
#'  the complete linkage clustering step of the RDP pipeline.
#'
#' @return An \code{\link{otuTable}} object parsed from the \code{".clust"} file. 
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
	return( otuTable(otumat, speciesAreRows=TRUE) )
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
#' @return An \code{otuTax} object containing both the otuTable and TaxonomyTable data
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
	# (otuTable and taxonomyTable)
	########################################
	### Initialize abundance matrix, a
	a <- matrix(0, nrow=length(x), ncol=(taxonomy_table_column_index-2))
	colnames(a) <- pyro_header[2:(taxonomy_table_column_index-1)]
	rownames(a) <- names(z)
	
	###### Initialize the raw pyrotagger taxonomy matrix, w
	ntaxtabcols <- (max(sapply(z, length)) - taxonomy_table_column_index + 1)
	w <- matrix("", nrow=length(x), ncol=ntaxtabcols)
	rownames(w) <- names(z)
	colnamesw <- pyro_header[-(1:(taxonomy_table_column_index-1))]
	colnamesw <- colnamesw[1:which(colnamesw=="Taxonomy")]
	colnamesw <- c(colnamesw, paste("col", (which(colnamesw=="Taxonomy")+1):ntaxtabcols, sep="") )
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
	OTU <- otuTable(a, speciesAreRows=TRUE)
	if( strict_taxonomy ){
		TAX <- taxTab[, biotaxonomy]
	} else {
		TAX <- taxTab(w)
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
#' Import mothur list and group files and return an otuTable
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
#' @return An \code{\link{otuTable}} object.
#'
#' @seealso \code{\link{import_mothur}}
#' @keywords internal
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
	# # original_tree <- read.tree(mothur_tree_file)
	# # tree          <- original_tree
	tree <- read.tree(mothur_tree_file)
	
	# Loop to merge the sequences in the tree by OTU
	# cycle through each otu, and sum the number of seqs observed for each sample
	for( i in names(otulist) ){
		# i <- names(otulist)[1]
		# First merge the reads that are in the same OTU ("eqspecies" argument)
		tree <- merge_species(tree, otulist[[i]])
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
#' @export
#' 
#' @examples
#' # # Load example data
#' # data(esophagus)
#' # export_env_file(esophagus, "~/Desktop/esophagus.txt")
export_env_file <- function(physeq, file="", writeTree=TRUE, return=FALSE){
	# data(esophagus)
	# physeq <- esophagus

	# Create otuTable matrix and force orientation
	OTU     <- as(otuTable(physeq), "matrix")
	if( !speciesAreRows(physeq) ){ OTU <- t(OTU) }
	
	# initialize sequence/sample names
	seqs    <- species.names(physeq)
	samples <- sample.names(physeq)	
	
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
		fileTree <- ps(file, ".nex")
		ape::write.nexus(tre(physeq), file=fileTree, original.data=FALSE)
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
# clip out the first 3 characters, and
# name according to the taxonomic rank
#' @keywords internal
parseGreenGenesPrefix <- function(char.vec){
	# Define the meaning of each prefix according to GreenGenes (and RDP?) taxonomy
	Tranks <- c(k="Kingdom", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species")
	taxvec        <- substr(char.vec, 4, 1000)
	names(taxvec) <- Tranks[substr(char.vec, 1, 1)]
	# Make sure order is same as Tranks
	taxvec <- taxvec[Tranks]; names(taxvec) <- Tranks
	return(taxvec)
}
################################################################################
#' Import phyloseq data from BIOM file
#'
#' New versions of QIIME produce a more-comprehensive and formally-defined
#' JSON file format. From the QIIME website:
#'
#' ``The biom file format (canonically pronounced `biome') is designed to be a 
#' general-use format for representing counts of observations in one or
#' more biological samples.''
#'
#' \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html} 
#'
#' @usage import_biom(BIOMfilename, taxaPrefix=NULL, parallel=FALSE, version=0.9)
#'
#' @param BIOMfilename (Required). A character string indicating the 
#'  file location of the BIOM formatted file. This is a JSON formatted file,
#'  specific to biological datasets, as described in 
#' 
#'  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#' 
#' @param taxaPrefix (Optional). Character string. What category of prefix precedes
#'  the taxonomic label at each taxonomic rank. Currently only ``greengenes'' is
#'  a supported option, and implies that the first letter indicates the 
#'  taxonomic rank, followed by two underscores and then the actual taxonomic
#'  assignment at that rank. The default value is \code{NULL}, meaning that
#'  no prefix or rank identifier will be interpreted. 
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
#' @return A \code{\link{phyloseq-class}} object.
#'
#' @seealso \code{\link{import}}, \code{\link{import_qiime}}
#'
#' @references \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#'
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr ldply
#' @importFrom plyr laply
#' @export
#' @examples
#'  # # # import with default parameters, specify a file
#'  # import_BIOM(myBIOMfile)
#'  # # # Example code for importing large file with parallel backend
#'  # library("doParallel")
#'  # registerDoParallel(cores=6)
#'  # import_biom("my/file/path/file.biom", taxaPrefix="greengenes", parallel=TRUE)
import_biom <- function(BIOMfilename, taxaPrefix=NULL, parallel=FALSE, version=0.9){
	
	# Read the data
	x <- fromJSON(BIOMfilename)
	
	########################################
	# OTU table:
	########################################
	# Check if sparse. Must parse differently than dense
	if( x$matrix_type == "sparse" ){
		otumat <- matrix(0, nrow=x$shape[1], ncol=x$shape[2])
		dummy <- sapply(x$data, function(i){otumat[(i[1]+1), (i[2]+1)] <<- i[3]})
	}
	# parse the dense matrix instead.
	if( x$matrix_type == "dense" ){
		# each row will be complete data values, should use laply
		# laply(x$data, vector, nrow=x$shape[1], ncol=x$shape[2])
		otumat <- laply(x$data, function(i){i}, .parallel=parallel)
	}
	
	# Get row (OTU) and col (sample) names
	rownames(otumat) <- sapply(x$rows, function(i){i$id})
	colnames(otumat) <- sapply(x$columns, function(i){i$id})
	
	otutab <- otuTable(otumat, TRUE)
	
	########################################
	# Taxonomy Table
	########################################
	# Need to check if taxonomy information is empty (minimal BIOM file)
	if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
		taxtab <- NULL
	} else {
		if( is.null(taxaPrefix) ){
			taxdf <- laply(x$rows, function(i){i$metadata$taxonomy}, .parallel=parallel)
		} else if( taxaPrefix == "greengenes" ){
			taxdf <- laply(x$rows, function(i){parseGreenGenesPrefix(i$metadata$taxonomy)}, .parallel=parallel)
		} else {
			taxdf <- laply(x$rows, function(i){i$metadata$taxonomy}, .parallel=parallel)
		}
		# Now convert to matrix, name the rows as "id" (the taxa name), coerce to taxonomyTable
		taxtab           <- as(taxdf, "matrix")
		rownames(taxtab) <- sapply(x$rows, function(i){i$id})
		taxtab <- taxTab(taxtab)	
	}
	
	########################################
	# Sample Data ("columns" in QIIME/BIOM)
	########################################
	# If there is no metadata (all NULL), then set samdata <- NULL
	if(  all( sapply(sapply(x$columns, function(i){i$metadata}), is.null) )  ){
		samdata <- NULL
	} else {
		samdata           <- ldply(x$columns, function(i){i$metadata}, .parallel=parallel)
		rownames(samdata) <- sapply(x$columns, function(i){i$id})
		samdata <- sampleData(samdata)
	}
	
	########################################
	# Put together into a phyloseq object
	########################################
	return( phyloseq(otutab, taxtab, samdata) )

}
################################################################################
################################################################################
################################################################################
################################################################################