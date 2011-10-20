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
#' @param envfilename A charater string of the ENV filename (relative or absolute)
#'
#' @param tree Optional phylo object that will be used to prune elements in the
#'  ENV hash table once it is read. Taxa not present as tip.labels in \code{tree}
#'  will be excluded from the output.
#'
#' @param sep A character string indicating the delimiter used in the file.
#'  The default is \code{"\t"}.
#'
#' @return A two- or three- column hash table (matrix).
#'
#' @keywords OTU trivial environment file
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
#' @param tipSampleTable a two-column character table (matrix or data.frame), 
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
# Another example.
# test2   <- read_env_file("~/Dropbox/ArrA/ML_full_190min_tradArrA_only_env_file.txt")
# testotu2<- envHash2otuTable( test2 )