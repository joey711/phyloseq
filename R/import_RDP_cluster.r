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