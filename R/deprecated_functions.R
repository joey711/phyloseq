# Deprecated function definitions and roxygen headers (source).
# Backward compatibility aliases



#' @rdname plot_richness
#' @aliases plot_richness
#' @export
plot_richness_estimates <- plot_richness


#' @rdname import_qiime_sample_data
#' @aliases import_qiime_sample_data
#' @export
import_qiime_sampleData <- import_qiime_sample_data

#' @rdname filterfun_sample
#' @aliases filterfun_sample
#' @export
filterfunSample <- filterfun_sample



#' @rdname genefilter_sample-methods
#' @aliases genefilter_sample
#' @export
genefilterSample <- genefilter_sample



#' @rdname prune_taxa-methods
#' @aliases prune_taxa
#' @export
prune_species <- prune_taxa



#' @rdname subset_taxa-methods
#' @docType methods
#' @export
subset_species <- subset_taxa



#' @rdname tip_glom-methods
#' @aliases tip_glom
#' @export
tipglom <- tip_glom



#' @rdname tax_glom
#' @aliases tax_glom
#' @export
taxglom <- tax_glom



#' @rdname phy_tree-methods
#' @aliases phy_tree
#' @export
tre <- phy_tree

#' @rdname assign-phy_tree
#' @aliases assign-phy_tree phy_tree<-
#' @export
"tre<-" <- function(x, value){
	phy_tree(x) <- value
	return(x)
}


#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows
#' @export
speciesarerows <- taxa_are_rows

#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows
#' @export
speciesAreRows <- taxa_are_rows

#' @rdname assign-taxa_are_rows
#' @aliases assign-taxa_are_rows
#' @export
"speciesAreRows<-" <- function(x, value){
	taxa_are_rows(x) <- value
	return(x)
}


#' @rdname ntaxa-methods
#' @aliases ntaxa
#' @export
nspecies <- ntaxa



#' @rdname taxa_names-methods
#' @aliases taxa_names
#' @export
species.names <- taxa_names



#' @rdname sample_names-methods
#' @aliases sample_names
#' @export
sampleNames <- sample_names

#' @rdname sample_names-methods
#' @aliases sample_names
#' @export
sample.names <- sample_names



#' @rdname get_sample-methods
#' @aliases get_sample
#' @export
getSamples <- get_sample



#' @rdname get_taxa-methods
#' @aliases get_taxa
#' @export
getSpecies <- get_taxa



#' @rdname rank_names
#' @aliases rank_names
#' @export
rank.names <- rank_names



#' @rdname get_taxa_unique
#' @aliases get_taxa_unique
#' @export
getTaxa <- get_taxa_unique



#' @rdname sample_variables
#' @aliases sample_variables
#' @export
sample.variables <- sample_variables



#' @rdname get_variable
#' @aliases get_variable
#' @export
getVariable <- get_variable



#' @rdname merge_taxa-methods
#' @aliases merge_taxa
#' @export
merge_species <- merge_taxa



#' @rdname otu_table-methods
#' @aliases otu_table
#' @export
otuTable <- otu_table

#' @rdname assign-otu_table
#' @aliases assign-otu_table otu_table<- otuTable<-
#' @export
"otuTable<-" <- function(x, value){
	otu_table(x) <- value
	return(x)
}



#' @rdname taxa_sums
#' @aliases taxa_sums
#' @export
speciesSums <- taxa_sums




#' @rdname sample_sums
#' @aliases sample_sums
#' @export
sampleSums <- sample_sums



#' @rdname sample_data-methods
#' @aliases sample_data
#' @export
sampleData <- sample_data

#' @rdname sample_data-methods
#' @aliases sample_data
#' @export
samData <- sample_data

#' @rdname assign-sample_data
#' @aliases sample_data<- 
#' @export
"sampleData<-" <- function(x, value){
	sample_data(x) <- value
	return(x)
}



#' @rdname tax_table-methods
#' @aliases tax_table
#' @export
taxtab <- tax_table

#' @rdname tax_table-methods
#' @aliases tax_table
#' @export
taxTab <- tax_table

#' @export
#' @rdname assign-tax_table
#' @aliases tax_table<-
"taxTab<-" <- function(x, value){
	tax_table(x) <- value
	return(x)
}



################################################################################
################################################################################
#' DEPRECATED. SEE \code{\link{psmelt}} converts OTU-table to data.frame
#'
#' Was used for plotting with \code{\link{plot_taxa_bar}} in the ggplot2 framework.
#'
#' @usage otu2df(physeq, taxavec, map, keepOnlyTheseTaxa=NULL, threshold=NULL)
#'
#' @param physeq An \code{otu_table} object.
#'
#' @param taxavec A character vector of the desired taxonomic names to 
#'  categorize each species in physeq. 
#' 
#' @param map The corresponding sample_data object for \code{physeq}.
#' 
#' @param keepOnlyTheseTaxa A vector of the taxonomic labels that you want
#'  included. If NULL, the default, then all taxonomic labels are used, except
#'  for the empty character string, ``'', which is trimmed away.
#'
#' @param threshold A [0,1] numeric. Fraction of abundance of the taxonomic
#'  groups to keep for each sample. The higher the value, the larger
#'  the diversity of taxonomic groups included. That is, a greater number of
#'  the rare groups are included. If NULL (or 1), the default, all taxonomic groups
#'  are included.
#'
#' @keywords internal
otu2df <- function(physeq, taxavec, map, keepOnlyTheseTaxa=NULL, threshold=NULL){
	########################################
	# sample_i - A sample name. A single character string.
	# physeq      - An otu_table object.
	# taxavec  - A character vector of the desired taxonomic names 
	#             to categorize each species in physeq. 
	otu2dfi <- function(sample_i, physeq, taxavec, normalized=TRUE){
		Abundance_i <- get_taxa(physeq, sample_i)
		if(normalized){ 
			Abundance_i <- Abundance_i / sum(Abundance_i)
		}
		dflist <- list(
			TaxaGroup = taxavec[taxa_names(physeq)],
			Abundance = Abundance_i,
			ID        = taxa_names(physeq),
			sample    = sample_i
		)
		dflist <- c(dflist, map[sample_i, , drop=TRUE])
		do.call("data.frame", dflist)
	}
	########################################
	# keepOnlyTheseTaxa - A character vector of the taxonomy labels to keep.
	# threshold     	- A [0, 1] numeric. Fraction of TaxaGroups to keep.
	#             the higher the value, the more TaxaGroups are included and
	#             the more the rare groups are included.
	# Trim subroutine
	trimdf <- function(df, keepOnlyTheseTaxa=NULL, threshold=NULL){
		# Force rows (OTUs) to be ordered by relative abundance
		df <- df[order(df$Abundance, decreasing=TRUE), ]
		# Trim nameless / unnecessary phyla from table:
		df <- subset(df, TaxaGroup != "", drop=TRUE)
		df <- subset(df, !is.na(TaxaGroup), drop=TRUE)
		# Keep a subset of phyla if keepOnlyTheseTaxa is specified
		if( !is.null(keepOnlyTheseTaxa) ){
			df <- subset(df, TaxaGroup %in% keepOnlyTheseTaxa, drop=TRUE)
		}
		# If threshold is provided, trim to most abundant threshold fraction.
		if( !is.null(threshold) ){
			top.TaxaGroup <- sort(tapply(df$Abundance, df$TaxaGroup,
								sum, na.rm=TRUE), decreasing=TRUE)
			top.TaxaGroup <- top.TaxaGroup[cumsum(top.TaxaGroup) <= threshold]
			df            <- subset(df, TaxaGroup %in% names(top.TaxaGroup))
		}
		# trim so that rows with zero abundance are removed
		df <- subset(df, Abundance > 0, drop=TRUE)
		return(df)		
	}
	########################################		
	# Main control loop to create large redundant data.frame for ggplot2	
	df <- otu2dfi(sample_names(physeq)[1], physeq, taxavec)
	df <- trimdf(df, keepOnlyTheseTaxa, threshold)
	for( j in sample_names(physeq)[-1] ){ 
		dfj <- otu2dfi( j, physeq, taxavec)
		dfj <- trimdf(dfj, keepOnlyTheseTaxa, threshold)
		df  <- rbind(df, dfj)
	}
	return(df)
}
################################################################################
#' DEPRECATED. USE \code{\link{plot_bar}} instead. Creates abundance barplot.
#'
#' This function wraps \code{ggplot2} plotting, and returns a \code{ggplot2}
#'  graphic object
#' that can be saved or further modified with additional layers, options, etc.
#' The main purpose of this function is to quickly and easily create informative
#' summary graphics of the differences in taxa abundance between samples in
#' an experiment. 
#'
#' The vertical axis is always relative abundance, but the data
#' can be further organized at the horizontal axis and faceting grid
#' by any combination of variates present in
#' the \code{\link{sample_data}} component of \code{physeq}.
#'
#' @usage plot_taxa_bar(physeq, taxavec="Domain",
#'	showOnlyTheseTaxa=NULL, threshold=NULL, x="sample", fill=x,  
#'	facet_formula = . ~ TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE, title=NULL)
#'
#' @param physeq (Required). An \code{\link{otu_table-class}} or 
#'  \code{\link{phyloseq-class}}.
#'  If \code{physeq} does not contain a taxonomyTable component,
#'  then the second argument, \code{taxavec}, is
#'  required and should have length equal to the number of species/taxa in
#'  \code{physeq}.
#'
#' @param taxavec (Optional, but recommended). A character vector of 
#'  the desired taxonomic rank to use in grouping 
#'  each species/OTU in \code{physeq}. If \code{physeq}  
#'  contains a taxonomyTable, then taxavec can alternatively specify the
#'  desired taxonomic level as a character string of length 1. 
#'  E.g. \code{taxavec = "Phylum"}. Default value is \code{"Domain"}.
#' 
#' @param showOnlyTheseTaxa A vector of the taxonomic labels that you want
#'  included. If NULL, the default, then all taxonomic labels are used, except
#'  for the empty character string, ``'', which is trimmed away.
#'
#' @param threshold (Optional). A [0,1] numeric. Fraction of abundance of the taxonomic
#'  groups to keep for each sample. The higher the value, the larger
#'  the diversity of taxonomica groups included. That is, a greater number of
#'  the rare groups are included. If NULL (or 1), the default, all taxonomic groups
#'  are included.
#'
#' @param x (Optional). A character string indicating which sample variable should be
#'  used to define the horizontal axis categories. Default is \code{"sample"}. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''.
#' 
#' @param fill (Optional). A character string. Indicates which sample variable
#'  should be used to define the fill color of the bars. This does not have to 
#'  match \code{x}, but does so by default. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''.
#' 
#' @param facet_formula A formula object as used by
#'  \code{\link{facet_grid}} in \code{\link{ggplot}} or \code{\link{qplot}}
#'  commands The default is: \code{. ~ TaxaGroup}. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''. E.g. An alternative
#'  \code{facet_grid} could be \code{sample ~ TaxaGroup}.
#'
#' @param OTUpoints (Optional). Logical. Default \code{FALSE}. Whether to add small grey 
#'  semi-transparent points for each OTU. Helps convey the relative distribution
#'  within each bar if it combines many different OTUs. For datasets with
#'  large numbers of samples and for complicated plotting arrangements, this
#'  might be too cluttered to be meaningful.
#'
#' @param labelOTUs (Optional). Logical. Default \code{FALSE}. Whether to add
#'  a label over the top
#'  few OTUs within each bar. As with \code{OTUpoints}, this is probably not
#'  a good idea for plots with large complexity. For low numbers of total OTUs
#'  this can be informative, and help display multiple layers of information 
#'  on the same graphic.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @return A ggplot2 graphic object.
#'
#' @seealso \code{\link{otu2df}}, \code{\link{qplot}}, \code{\link{ggplot}}
#'
#' @import ggplot2
#' @export
#' @aliases taxaplot
#' @rdname plot-taxa-bar
#'
#' @examples
#' data("GlobalPatterns")
#' gp.ch = subset_species(GlobalPatterns, Phylum == "Chlamydiae")
#' plot_taxa_bar(gp.ch, "Genus", x="SampleType", fill="SampleType", title="plot_taxa_bar: deprecated. Do not use.")
plot_taxa_bar <- function(physeq, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x="sample", fill=x, 
	facet_formula = . ~ TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE, title=NULL){

	# Some preliminary assignments. Assumes physeq has non-empty sample_data slot.
	map <- sample_data(physeq)
	if( length(taxavec) == 1 ){ 
		taxavec <- as(tax_table(physeq), "matrix")[, taxavec, drop=TRUE]
	}

	# Build the main species-level data.frame
	df <- otu2df(physeq, taxavec, map,	showOnlyTheseTaxa, threshold)

	########################################
	# Set the factor-order for df to ensure the
	# taxonomic groups are neatly ordered
	# by their relative contribution in the total dataset.
	top.TaxaGroup   <- sort(tapply(df$Abundance, df$TaxaGroup,
						sum, na.rm=TRUE), decreasing=TRUE)
	top.TaxaGroup   <- 	top.TaxaGroup / sum(top.TaxaGroup)
	top.TaxaGroup   <- 	top.TaxaGroup[top.TaxaGroup > 0]						
	   df$TaxaGroup <- factor(as.character(df$TaxaGroup),
						levels=names(top.TaxaGroup))

	# If there is a showOnlyTheseTaxa provided, use that as the order instead.
	if( !is.null(showOnlyTheseTaxa) ){
		df$TaxaGroup <- factor(df$TaxaGroup, levels=showOnlyTheseTaxa)	
	}

	########################################
	# Create taxaGroup-summed df for making solid bars, dftot
	df2sampleTGtot <- function(df, map){
		AbundTot <- tapply(df$Abundance, list(df$sample, df$TaxaGroup), sum)
		
		TaxaGrouptot <- as.vector(sapply( levels(df$TaxaGroup),
						rep, times = length(levels(df$sample)) ))
		sampletot    <- rep(levels(df$sample), times=length(levels(df$TaxaGroup)))
		
		dftot <- data.frame(
			TaxaGroup = TaxaGrouptot,
			sample    = sampletot,
			Abundance = as.vector(AbundTot)
		)
		dftot <- cbind(dftot, map[factor(dftot$sample), ])
		# Make sure TaxaGroup levels-order matches df
		dftot$TaxaGroup <- factor(dftot$TaxaGroup, levels(df$TaxaGroup))		
		return(dftot)
	}
	dftot <- df2sampleTGtot(df, map)

	########################################
	# Build the ggplot
	p  <- ggplot(df) + 
		theme(axis.text.x=element_text(angle=-90, hjust=0))

	p <- p + 
		# The full stack
		geom_bar(
			data=dftot, 
			eval(call("aes",
				x=as.name(x), 
				y=quote(Abundance), 
				fill=as.name(fill)
			)),
			position="dodge", stat="identity"
		) + 
		# Some reasonable default options
		theme(panel.grid.minor = element_blank()) + 
		theme(panel.grid.major = element_blank()) +
		theme(panel.border = element_blank()) +
		labs(y="Relative Abundance", x=x, fill=fill)
		
	# Should the individual OTU points be added. Default FALSE
	if( OTUpoints ){
		p <- p + geom_point(
			data=df, 
			eval(call("aes",
				x=as.name(x), 
				y=quote(Abundance)
			)),
			color="black", size=1.5, position="jitter", alpha=I(1/2)
		)
	}
		
	# Should the most abundant OTUs be labeled. Default FALSE
	if( labelOTUs ){
		# Create a small df subset for labelling abundant OTUs
		dfLabel <- subset(df, Abundance > 0.05)	
		p <- p + geom_text(data=dfLabel, size=2,
			eval(call("aes", 
				x=as.name(x), 
				y=quote(Abundance+0.01), 
				label=quote(ID)
			)),
		)
	}
		

	if( !is.null(facet_formula) ){	
		p <- p + facet_grid(facet_formula)
	}
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}	
	
	########################################
	# Return the ggplot object so the user can 
	# additionally manipulate it.
	return(p)
}
################################################################################
#' @export
#' @aliases plot_taxa_bar
#' @rdname plot-taxa-bar
taxaplot <- plot_taxa_bar
################################################################################
################################################################################


