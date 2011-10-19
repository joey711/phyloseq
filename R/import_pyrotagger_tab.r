################################################################################
#' Imports a tab-delimited version of the pyrotagger output file.
#'
#' PyroTagger is a web-server that takes raw, barcoded 16S rRNA amplicon sequences
#' and returns an excel spreadsheet (\code{".xls"}) with both abundance and 
#' taxonomy data. It also includes some confidence information related to the 
#' taxonomic assignment.
#'
#' PyroTagger is served by the Joint Genome Institute at \code{"http://pyrotagger.jgi-psf.org/"}
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