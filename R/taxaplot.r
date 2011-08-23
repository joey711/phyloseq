################################################################################
#' Convert an otuTable object into a redundant data.frame useful for plotting
#' in the ggplot2 framework.
#'
#' @param otu An \code{otuTable} object.
#'
#' @param taxavec A character vector of the desired taxonomic names to 
#'  categorize each species in otu. 
#' 
#' @param map The corresponding sampleMap object for \code{otu}.
#' 
#' @param keepOnlyTheseTaxa A vector of the taxonomic labels that you want
#'  included. If NULL, the default, then all taxonomic labels are used, except
#'  for the empty character string, ``'', which is trimmed away.
#'
#' @param threshold A [0,1] numeric. Fraction of abundance of the taxonomic
#'  groups to keep for each sample. The higher the value, the larger
#'  the diversity of taxonomica groups included. That is, a greater number of
#'  the rare groups are included. If NULL (or 1), the default, all taxonomic groups
#'  are included.
#'
#' @seealso taxaplot
#' @export
#' @examples #
otu2df <- function(otu, taxavec, map, keepOnlyTheseTaxa=NULL, threshold=NULL){
	########################################
	# sample_i - A sample name. A single character string.
	# otu      - An otuTable object.
	# taxavec  - A character vector of the desired taxonomic names 
	#             to categorize each species in otu. 
	otu2dfi <- function(sample_i, otu, taxavec, normalized=TRUE){
		Abundance_i <- getSpecies(otu, sample_i)
		if(normalized){ 
			Abundance_i <- Abundance_i / sum(Abundance_i)
		}
		dflist <- list(
			TaxaGroup = taxavec[species.names(otu)],
			Abundance = Abundance_i,
			ID        = species.names(otu),
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
	df <- otu2dfi(sample.names(otu)[1], otu, taxavec)
	df <- trimdf(df, keepOnlyTheseTaxa, threshold)
	for( j in sample.names(otu)[-1] ){ 
		dfj <- otu2dfi( j, otu, taxavec)
		dfj <- trimdf(dfj, keepOnlyTheseTaxa, threshold)
		df  <- rbind(df, dfj)
	}
	return(df)
}
################################################################################
#' Create a structured barplot graphic of the taxonomic groups.
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
#' the sampleMap component of \code{otu}, or the \code{map} argument if
#' \code{otu} is a simple \code{otuTable}.
#'
#' @param otu An \code{otuTable} object, or higher-order object that contains
#'  an otuTable and sampleMap (e.g. ``phyloseq'' class and its superclasses.)
#'
#' @param taxavec A character vector of the desired taxonomic names to 
#'  categorize each species in \code{otu}. If \code{otu} is a higher-order
#'  object that
#'  contains a taxonomyTable, then taxavec can alternatively specify the
#'  desired taxonomic level as a character string of length 1. 
#'  E.g. \code{taxavec = "Phylum"}. Default value is \code{"Domain"}.
#' 
#' @param map The corresponding sampleMap object for \code{otu}. This is ignored
#'  for higher-order objects, where the sampleMap component is used instead.
#' 
#' @param showOnlyTheseTaxa A vector of the taxonomic labels that you want
#'  included. If NULL, the default, then all taxonomic labels are used, except
#'  for the empty character string, ``'', which is trimmed away.
#'
#' @param threshold A [0,1] numeric. Fraction of abundance of the taxonomic
#'  groups to keep for each sample. The higher the value, the larger
#'  the diversity of taxonomica groups included. That is, a greater number of
#'  the rare groups are included. If NULL (or 1), the default, all taxonomic groups
#'  are included.
#'
#' @param \code{x_category} A character string indicating which map column should be
#'  used to define the horizontal axis categories. Default is \code{"sample"}. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''.
#' 
#' @param \code{fill_category} A character string indicating which map column
#'  should be used to define the fill color of the bars. This does not have to 
#'  match \code{x_category}, but does so by default. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''.
#' 
#' @param \code{facet_formula} A formula object as used by
#'  \code{\link{facet_grid}} in \code{\link{ggplot}} or \code{\link{qplot}}
#'  commands The default is: \code{. ~ TaxaGroup}. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''. E.g. An alternative
#'  \code{facet_grid} could be \code{sample ~ TaxaGroup}.
#'
#' @seealso \code{\link{otu2df}}, \code{\link{qplot}}, \code{\link{ggplot}}
#' @export
#' @examples #
#' # data(ex1)
#' # taxaplot(ex1, "Class", threshold=0.85, x_category="Diet",
#' # fill_category="Diet", facet_formula = Gender ~ TaxaGroup)
setGeneric("taxaplot", function(otu, taxavec="Domain", map,
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category,  
	facet_formula = . ~ TaxaGroup) standardGeneric("taxaplot"))
# The main method.
setMethod("taxaplot", "otuTable", function(otu, taxavec="Domain", map,
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category, 
	facet_formula = . ~ TaxaGroup){

	# Build the main species-level data.frame
	df <- otu2df(otu, taxavec, map,	showOnlyTheseTaxa, threshold)

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

	# Create a small df subset for labelling abundant OTUs
	dfLabel <- subset(df, Abundance > 0.05)	

	########################################
	# Build the ggplot
	p  <- ggplot(df) + 
		opts(axis.text.x=theme_text(angle=-90, hjust=0))

	p <- p + 
			# geom_bar(
			# data=df, eval(call("aes", x=as.name(x_category), y=quote(Abundance), fill=as.name(fill_category) )),
			#	 position="stack", stat="identity"
			# ) +
		# The full stack
		geom_bar(
			data=dftot, 
			eval(call("aes",
				x=as.name(x_category), 
				y=quote(Abundance), 
				fill=as.name(fill_category)
			)),
			position="dodge", stat="identity"
		) + 
		geom_point(
			data=df, 
			eval(call("aes",
				x=as.name(x_category), 
				y=quote(Abundance)
			)),
			color="black", size=1.5, position="jitter", alpha=I(1/2)
		) +
		opts(panel.grid.minor = theme_blank()) + 
		opts(panel.grid.major = theme_blank()) +
		opts(panel.border = theme_blank()) +
		geom_text(
			data=dfLabel,
			size=2,
			eval(call("aes", 
				x=as.name(x_category), 
				y=quote(Abundance+0.01), 
				label=quote(ID)
			)),
		) +	
		labs(y="Relative Abundance", x=x_category, fill=fill_category)

	if( !is.null(facet_formula) ){	
		p <- p + facet_grid(facet_formula)
	}
	########################################
	# Return the ggplot object so the user can 
	# additionally manipulate it.
	return(p)
# Close otuTable taxaplot method.
})
setMethod("taxaplot", "phyloseq", function(otu, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category, 
	facet_formula= . ~ TaxaGroup){

	taxaplot(otuTable(otu), taxavec, map=sampleMap(otu),
		showOnlyTheseTaxa, threshold,  x_category, fill_category, facet_formula)
})
setMethod("taxaplot", "phyloseqTax", function(otu, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category, 
	facet_formula= . ~ TaxaGroup){

	if( length(taxavec) == 1 ){ 
		taxavec <- as(taxTab(otu), "matrix")[, taxavec, drop=TRUE]
	}

	taxaplot(phyloseq(otu), taxavec=taxavec, showOnlyTheseTaxa=showOnlyTheseTaxa,
		threshold=threshold, x_category=x_category, fill_category=fill_category,
		facet_formula=facet_formula)
})
################################################################################