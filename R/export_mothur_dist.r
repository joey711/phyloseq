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
#' ### data(ex1) 
#' ### myDistObject <- as.dist(cophenetic(tre(ex1)))
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