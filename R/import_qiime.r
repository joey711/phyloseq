######################################################################################
#' Import function to read files created by the QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' including especially an OTU file that typically contains both OTU-abundance
#' and taxonomic identity information. The map-file is also an important input
#' to QIIME that stores sample covariates, converted naturally to the 
#' sampleMap component data type in the phyloseq-package. QIIME may also produce a
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
#'  ()\code{"phylo4"} or \code{"phylo"} class), in case the tree has already been
#'  imported, or is in a different format than NEXUS.
#'
#' @param biotaxonomy (Optional). A character vector indicating the name of each taxonomic level
#'  in the taxonomy-portion of the otu-file, which may not specify these levels 
#'  explicitly.
#'  Default value is \code{NULL}. 
#'
#' @param ... Additional arguments passed to \code{\link{read.nexus}}, as necessary.
#'  Make sure that your phylogenetic tree file is readable by \code{\link{read.nexus}}
#'  prior to calling this function.
#'
#' @return A \code{\link{phyloseq-class}} object.
#'
#' @seealso \code{\link{phyloseq}}, \code{\link{merge_phyloseq}},
#'	 \code{\link{read.tree}}, \code{\link{read.nexus}}, \code{\link{readNexus}}
#'   \code{\link{readNewick}}
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
#' @examples #
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
		QiimeMap     <- import_qiime_sampleMap(mapfilename)
		argumentlist <- c(argumentlist, list(QiimeMap))
	}

	if( !is.null(otufilename) ){
		otutax       <- import_qiime_otu_tax(otufilename, biotaxonomy)	
		argumentlist <- c(argumentlist, list(otuTable(otutax)), list(taxTab(otutax)) )
	}

	if( !is.null(treefilename) ){
		if( class(treefilename) %in% c("phylo", "phylo4") ){ 
			tree <- treefilename
		} else {
			#tree <- read.tree(treefilename, ...)
			#tree <- readNexus(treefilename, ...)
			#tree <- readNewick(treefilename, ...)
			tree <- ape::read.nexus(treefilename, ...)
		}
		if( class(tree)=="phylo"){
			tree <- as(tree, "phylo4")
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
#'   \code{\link{import_qiime_sampleMap}}
#'
#' @keywords internal
#' @examples #
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
#' Import just \code{sampleMap} file from QIIME pipeline.
#'
#' QIIME produces several files that can be analyzed in the phyloseq-package, 
#' This includes the map-file, which is an important \emph{input}
#' to QIIME that can also indicate sample covariates. It is converted naturally to the 
#' sampleMap component data type in phyloseq-package, based on the R data.frame.   
#' 
#' See \code{\link{import_qiime}} for more information about QIIME. It is also the
#' suggested function for importing QIIME-produced data files. 
#'
#' @usage import_qiime_sampleMap(mapfilename)
#'
#' @param mapfilename (Required). A character string. The name of the QIIME map
#'  file required for processing pyrosequencing tags
#'  in QIIME as well as some of the post-clustering analysis. This is a required
#'  input file for running QIIME. Its strict formatting specification is expected by
#'  this function, do not attempt to modify it manually once it has worked properly
#'  in QIIME. 
#'
#' @return A \code{sampleMap} object.
#'
#' @seealso \code{\link{import_qiime}}, \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
#'	 \code{\link{import_qiime_otu_tax}}
#'
#' @keywords internal
#' @examples #
import_qiime_sampleMap <- function(mapfilename){
	# Process mapfile. Name rows as samples.
	QiimeMap <- read.table(file=mapfilename, header=TRUE,
		sep="\t", comment.char="")
	rownames(QiimeMap) <- as.character(QiimeMap[,1])
	return( sampleMap(QiimeMap) )
}
######################################################################################
################################################################################
#' Read a UniFrac-formatted ENV file.
#'
#' Convenience wrapper function to read the environment-file, as formatted for 
#' input to the UniFrac server(REF). The official format of these files is that
#' each row specifies (in order) the sequence name, source sample, and (optionally)
#' the number of times the sequence was observed. 
#'
#' @usage read_env_file(envfilename, tree=NULL, sep="\t")
#'
#' @param envfilename (Required). A charater string of the ENV filename (relative or absolute)
#'
#' @param tree (Optional). A phylo object that will be used to prune elements in the
#'  ENV hash table once it is read. Taxa not present as tip.labels in \code{tree}
#'  will be excluded from the output.
#'
#' @param sep A character string indicating the delimiter used in the file.
#'  The default is \code{"\t"}.
#'
#' @return A two- or three- column hash table (matrix).
#'
#' @export
#' @examples #
read_env_file <- function(envfilename, tree=NULL, sep="\t"){
	tipSampleTable <- read.table(envfilename, sep=sep)
	# Remove elements of tipSampleTable that aren't in the tree
	if( class(tree) == "phylo" ){
		tipSampleTable <- tipSampleTable[tipSampleTable[,1] %in% tree$tip.label,]
	}
	return(tipSampleTable)
}
################################################################################  
#' Convert a sequence-sample hash (like ENV file) into an OTU table.
#' 
#' Takes a two- or three-column character table (matrix or data.frame), where each
#' row specifies the sequence name, source sample, and (optionally) abundance. 
#' This format mirrors the output
#' from \code{\link{read_env_file}}, consistent with the env file for UniFrac
#' server. Parses the 2 column table into a sparse matrix of species-by-sample, where
#' each species-row has only one non-zero value. We call this sparse abundance
#' table the trivial OTU table, where every sequence is treated as a separate 
#' species. If a phylogenetic tree is available, it can be submitted with this
#' table as arguments to \code{\link{tipglom}} to create an object with a
#' non-trivial \code{otuTable}.  
#'
#' @usage envHash2otuTable(tipSampleTable)
#'
#' @param tipSampleTable (Required). A two-column character table (matrix or data.frame), 
#' where each row specifies the sequence name and source sample. This format 
#' mirrors the output from \code{\link{read_env_file}}, consistent with the 
#' env-file for the UniFrac server(REF). 
#'
#' @return Trivial (sparse) OTU table. Object of class \code{otuTable}. 
#' The 2- (or 3-) column hash table is parsed into a sparse matrix of 
#' species-by-sample, where
#' each species-row has only one non-zero value. We call this sparse abundance
#' table the trivial OTU table, where every sequence is treated as a separate 
#' species.
#' 
#' @seealso tipglom read_env_file otuTable
#'
#' @keywords OTU trivial environment file
#' @export
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
	otuTable(trivialOTU, speciesAreRows=TRUE)
}
envhash2otutable <- envHash2otuTable
################################################################################