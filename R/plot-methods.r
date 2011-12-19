#
# extension of plot methods for phyloseq object.
# 
################################################################################
################################################################################
################################################################################
################################################################################
#' Generic plot defaults for phyloseq.
#'
#' The specific plot type is chosen according to available non-empty slots.
#' This is mainly for syntactic convenience and quick-plotting. See links below
#' for some examples of available graphics tools available in the
#' \code{\link{phyloseq-package}}.  
#'
#' @usage plot_phyloseq(physeq, ...)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}. The actual plot type
#'  depends on the available (non-empty) component data types contained within.
#'
#' @param ... (Optional). Additional parameters to be passed on to the respective
#'  specific plotting function. See below for different plotting functions that
#'  might be called by this generic plotting wrapper.
#'
#' @return A plot is created. The nature and class of the plot depends on 
#'  the \code{physeq} argument, specifically, which component data classes
#'  are present.
#' 
#' @seealso 
#'  \code{\link{taxaplot}}
#'  \code{\link{plot_tree_phyloseq}}
#'  \code{\link{plot_ordination_phyloseq}}
#'  \code{\link{calcplot}}
#'  \code{\link{makenetwork}}
#'
#' @export
#' @docType methods
#' @rdname plot_phyloseq-methods
#'
#' @examples 
#'  ## data(ex1)
#'  ## plot_phyloseq(ex1)
setGeneric("plot_phyloseq", function(physeq, ...){ standardGeneric("plot_phyloseq") })
#' @aliases plot_phyloseq,phyloseq-method
#' @rdname plot_phyloseq-methods
setMethod("plot_phyloseq", "phyloseq", function(physeq, ...){
	if( all(c("otuTable", "sampleData", "tre") %in% getslots.phyloseq(physeq)) ){
		plot_tree_phyloseq(physeq, ...)		
	} else if( all(c("otuTable", "sampleData", "taxTab") %in% getslots.phyloseq(physeq) ) ){
		taxaplot(otu=physeq, ...)
	} else if( all(c("otuTable", "tre") %in% getslots.phyloseq(physeq)) ){
		tree <- tre(physeq)
		ape::plot.phylo(tree, ...)	
		ape::nodelabels(as.character(1:max(tree$edge)), node=1:max(tree$edge))
		ape::edgelabels(as.character(1:nrow(tree$edge)), edge=1:nrow(tree$edge))		
	} else {
		makenetwork(physeq)
	}
})
################################################################################
################################################################################
# calcplotrda and calcplotcca
################################################################################
#' Convenient rendering of ordination results using ggplot2.
#'
#' This convenience function includes many useful defaults for generating a
#' plot of your ordination. This supplements the default plots created by the
#' plot method extension from the vegan package, \code{\link{plot.cca}}.
#' For further examples, see the phyloseq vignette.
#'
#' Because the ggplot2 objects can be stored as variables, and then further 
#' modified, it is not necessary to support additional arguments as \code{...}.
#' If further modifications to the plot produced by this function are desired,
#' they can be provided as parameter updates to the returned object. See 
#' the specific examples provided by the phyloseq vignette, or the more 
#' comprehensive documentation available for the ggplot2 package itself. 
#' 
#' @usage plot_ordination_phyloseq(mod, object,
#'	plot_title = as(mod$call, "character")[1],
#'	man.colors=NULL, man.sizes=NULL, man.shapes=NULL,
#'	species_alpha=1/10, sites_alpha=2/3,
#'	species_color_category =NULL,
#'	species_shape_category =NULL,
#'	species_size_category  =NULL,
#'	site_color_category=NULL,
#'	site_shape_category=NULL,
#'	site_size_category=NULL,	
#'	add_sites_data = NULL,
#'	add_taxa_data = NULL)
#'
#' @param mod (Required). A \code{cca} or \code{rda} results object. 
#'  See \code{\link{cca.object}}. For phyloseq
#'  objects, you probably want to see \code{\link{cca.phyloseq}} or 
#'  \code{\link{rda.phyloseq}}.
#'  
#' @param object (required). A \code{\link{phyloseq-class}} object. Necessary
#'  because ordination result objects
#'  will only refer to the abundance table, not the phyloseq object used to
#'  create them.
#' 
#' @param plot_title Default is to use the name of the ordination method, stored in 
#'  the ordination result object, accessed via \code{as(mod$call, "character")[1]}.
#'  Any character string will suffice.
#' 
#' @param man.colors A character vector of colors to overide the default 
#'  color scale used by ggplot2. Default is \code{NULL}. 
#' 
#' @param man.sizes An integer vector of sizes to overide the default size
#'  scale used by ggplot2. Default is \code{NULL}. 
#' 
#' @param man.shapes An integer vector of shape indices to overide the defaults
#'  used by ggplot2. Default is \code{NULL}. 
#' 
#' @param species_alpha A single numeric. The alpha (transparency) value
#'  to use for all \code{geom_points}
#'  that represent taxa. Default is \code{1/10}.
#' 
#' @param sites_alpha A single numeric. The alpha (transparency) value
#'  to use for all \code{geom_points}
#'  that represent samples/sites. Default is \code{2/3}.
#' 
#' @param species_color_category A manual color scale for taxa.
#'  Default is \code{NULL}. 
#'
#' @param species_shape_category A manual shape scale for taxa.
#' Default is \code{NULL}. 
#'
#' @param species_size_category A manual size scale for taxa.
#' Default is \code{NULL}.
#' 
#' @param site_color_category A manual color scale for samples/sites.
#'  Default is \code{NULL}. 
#'
#' @param site_shape_category A manual shape scale for samples/sites.
#'  Default is \code{NULL}.
#'
#' @param site_size_category A manual size scale for samples/sites.
#'  Default is \code{NULL}.
#'
#' @param add_sites_data A \code{data.frame} object providing additional data
#'  about the samples/sites that is not already available in the \code{sampleData} of
#'  your phyloseq object specified by \code{object}. As such, it must have 
#'  the same number of rows as as there are samples in \code{object}. 
#'  Default is \code{NULL}. 
#' 
#' @param add_taxa_data A \code{matrix} object providing additional data
#'  about the taxa that is not already available in the \code{taxTab} of
#'  your phyloseq object specified by \code{object}, if it has a \code{taxTab} slot.
#'  It is not required. As such, it must have 
#'  the same number of rows as there are taxa in \code{object}. 
#'  Default is \code{NULL}.
#'
#' @return A graphic object from the ggplot2 package.
#'
#' @seealso \code{\link{calcplot}}, \code{\link[vegan]{plot.cca}}, \code{\link[vegan]{cca.object}}
#' @import vegan
#' @export
#' @examples
#' # data(ex1)
#' # ex4  <- transformsamplecounts(ex1, threshrankfun(500))
#' # # RDA
#' # modr <- rda.phyloseq(ex4 ~ Diet + Gender)
#' # # CCA
#' # modc <- cca.phyloseq(ex1 ~ Diet + Gender)
#' # plot_ordination_phyloseq(modr, ex1)
#' # plot_ordination_phyloseq(modc, ex1)
#' # plot_ordination_phyloseq(modr, ex1, species_color_category="Phylum", species_alpha=1/5)
#' # plot_ordination_phyloseq(modr, ex1, species_color_category="Phylum",
#' # site_shape_category="Diet", site_color_category="Gender", species_alpha=1)
#' # plot_ordination_phyloseq(modr, ex1, species_color_category="Phylum",
#' # site_shape_category="Diet", site_color_category="Gender", site_size_category="total.reads",
#' # species_alpha=0.4, sites_alpha=2/3,
#' # add_sites_data=data.frame(total.reads=sampleSums(ex1)) 
#' # site.colors <- c("darkorchid1", "blue3")
#' # species.col.groups <- unique(as(taxTab(ex1), "matrix")[,"Phylum"])
#' # man.colors <- c(site.colors, rainbow((1+length(species.col.groups)),start=11/12, end=5/12))
#' # p<-plot_ordination_phyloseq(modr, ex1, species_color_category="Phylum",
#' # site_shape_category="Diet", site_color_category="Gender", site_size_category="total.reads",
#' # species_alpha=0.4, sites_alpha=2/3,
#' # add_sites_data=data.frame(total.reads=sampleSums(ex1)),
#' # man.colors=man.colors)
#' # print(p)
plot_ordination_phyloseq <- function(mod, object,
	plot_title = as(mod$call, "character")[1],
	man.colors=NULL, man.sizes=NULL, man.shapes=NULL,
	species_alpha=1/10, sites_alpha=2/3,
	species_color_category =NULL,
	species_shape_category =NULL,
	species_size_category  =NULL,
	site_color_category=NULL,
	site_shape_category=NULL,
	site_size_category=NULL,	
	add_sites_data = NULL,
	add_taxa_data = NULL
			){
	#require("ggplot2")

	########################################
	# Create sites and species (taxa) data.frames
	########################################
	pmod      <- vegan::scores(mod, display=c("biplot","cn","sites","species"), scaling=1)
	sitesdata <- data.frame(axis1 = pmod$sites[, 1],   axis2 = pmod$sites[, 2]   )
	specidata <- data.frame(axis1 = pmod$species[, 1], axis2 = pmod$species[, 2] )

	########################################
	# If there is additional data in the otuSam-class
	# object specified as argument 'object' add it
	# automatically to the respective data.frame
	########################################
	# if there is a sampleData, cbind it to the sites data
	if( !is.null(sampleData(object)) ){
		sitesdata <- cbind(sitesdata, data.frame(sampleData(object)))
	}
	
	# if there is a taxTab, cbind it to the species data
	if( !is.null(taxTab(object)) ){
		specidata <- cbind(specidata, data.frame(taxTab(object)))
	}

	########################################
	# if have additional species or sites data to include, cbind-it
	########################################
	if( !is.null(add_sites_data) ){
		sitesdata <- cbind(sitesdata, add_sites_data)
	}
	if( !is.null(add_taxa_data) ){
		specidata <- cbind(specidata, add_taxa_data)
	}
	
	########################################
	# setup ggplot sitesdata first.
	########################################
	# First, build the argument list.
	arglist=list(x=quote(axis1), y=quote(axis2))
	if( !is.null(site_color_category) ){
		arglist <- c(arglist, list(colour=as.name(site_color_category)))
	}
	if( !is.null(site_shape_category) ){
		arglist <- c(arglist, list(shape =as.name(site_shape_category)))
	}
	if( !is.null(site_size_category) ){
		arglist <- c(arglist, list(size  =as.name(site_size_category)))
	}	
	# save sites_arglist for re-plot at end
	sites_arglist <- arglist

	# Finally, initalize ggplot object as "p"
	p <- ggplot2::ggplot(data=sitesdata)
	# Start plot with dummy layer, invisible sites
	p <- p + ggplot2::geom_point( do.call("aes", sites_arglist), alpha=0)

	########################################
	# Add manual values for color, size, and shape.
	########################################
	if( !is.null(man.colors) ){
		p <- p + ggplot2::scale_color_manual(value = man.colors)
	}
	if( !is.null(man.shapes) ){
		p <- p + ggplot2::scale_shape_manual(values = man.shapes)
	}
	if( !is.null(man.sizes) ){
		p <- p + ggplot2::scale_size_manual(values = man.sizes)
	}
	
	########################################
	## Add the species data layer. 
	########################################
	# First, make the argument list for site aesthetic
	arglist=list(x=quote(axis1), y=quote(axis2))
	if( !is.null(species_color_category) ){
		arglist <- c(arglist, list(colour=as.name(species_color_category)))
	}
	if( !is.null(species_shape_category) ){
		arglist <- c(arglist, list(shape =as.name(species_shape_category)))
	}	
	if( !is.null(species_size_category) ){
		arglist <- c(arglist, list(size  =as.name(species_size_category)))
	}
	p <- p + ggplot2::geom_point(data=specidata, do.call("aes", arglist), alpha=species_alpha)

	########################################
	## Re-Add the sites-layer after the species, 
	## since there will normally be many fewer of them than species.
	########################################
	if( is.null(site_size_category) ){
		p <- p + ggplot2::geom_point( do.call("aes", sites_arglist), size=7, alpha=sites_alpha)		
	} else {
		p <- p + ggplot2::geom_point( do.call("aes", sites_arglist), alpha=sites_alpha)		
	}

	########################################
	# Adjust the background grey slightly so that it is more
	# easily visible on electronic displays
	########################################
	p <- p + ggplot2::opts(panel.background = ggplot2::theme_rect(fill = "grey76", colour = NA),
		strip.background = ggplot2::theme_rect(fill = "grey76", colour = NA)
	)	
	
	########################################
	# Add a plot title based on the ordination type in mod
	########################################
	p <- p + ggplot2::opts(title = plot_title)
	
	########################################
	# Adjust the legend title for color.
	# Mixed legends default to the first title,
	# which won't make sense. Rename it.
	########################################
	if( !is.null(species_color_category) & !is.null(site_color_category) ){
		p <- p + ggplot2::labs(colour = "Sites, Species Colors")
	}
	
	return(p)
}
################################################################################
# calcplot
################################################################################
#' Convenience wrapper for performing ordination and plotting.
#'
#' @usage calcplot(X, RDA_or_CCA="cca", object=get(all.vars(X)[1]), ...)
#'
#' @param X (Required). A formula object. The left-hand side specifying a single 
#'  phyloseq object that contains (at minimum) an \code{otuTable} and a 
#'  \code{sampleData}. The right-hand side should contain the label of at least
#'  one variate present in the \code{sampleData} of the LHS object. There is
#'  one supported alternative, the special \code{"."} for the RHS, in which case
#'  all available variates will be specified as constraints in the constrained
#'  ordination, but not shaded. This behavior is less explicit and may depend
#'  somewhat arbitrarily on the order of variates in your \code{sampelMap}; 
#'  as such, it should be used with caution.
#'
#'  For example, there are two variates in the example object obtained with
#'  \code{data(ex1)}. \code{X} can be specified as 
#'
#'  \code{ex1 ~ Diet + Gender}
#'
#'  for using both as constraints in the ordination (and shading), or as
#'
#'  \code{ex1 ~ Diet}
#'
#'  if you only wanted to constrain on the Diet variable, for example.
#'
#'  For available
#'  variate names in your object, try \code{colnames(sampleData(ex1))}, where
#'  \code{ex1} is the phyloseq object containing your data. Because this is
#'  a formula object, quotes should not be used. See \code{\link{formula}}
#'  for details about writing a formula in \code{R}.
#'
#' @param RDA_or_CCA A character string, indicating whether the ordination
#'  method should be RDA or CCA. Default is "cca". Case ignored. Only first
#'  letter of string is considered.
#'
#' @param object Default is the object specified in the left-hand 
#'  side of \code{X}.
#'
#' @param ... Additional plotting arguments, passed on to 
#'  \code{\link{plot_ordination_phyloseq}}
#'
#' @return A ggplot2 graphics object. If not stored as a variable, a graphic
#'  object will be produced on the default device. 
#'
#' @export
#' @seealso \code{\link{rda.phyloseq}}, \code{\link{cca.phyloseq}},
#'  \code{\link{rda}}, \code{\link{cca}}
#' @examples
#' ## data(ex1)
#' ## calcplot(ex1 ~ Diet + Gender)
calcplot <- function(X, RDA_or_CCA="cca", object=get(all.vars(X)[1]), ...){
	if( class(X) != "formula" ){
		warning("First argument, X, must be formula class. See ?formula for more info\n")
		return()
	}
	
	if( substr(RDA_or_CCA, 1, 1) %in% c("R", "r") ){
		mod <- rda.phyloseq(X)		
	} else if( substr(RDA_or_CCA, 1, 1) %in% c("C", "c") ){
		mod <- cca.phyloseq(X)		
	} else {
		cat("You did not properly specify the desired ordination method\n")
		cat("Please see documentation.\n")		
		return()
	}

	ord_vars <- all.vars(X)[-1]

	# Initialize call list for plot_ordination_phyloseq
	popcallList <- list(mod=mod, object=object)

	# Populate the shading/shaping options in the call list.
	if( ord_vars[1] != "."){
		popcallList <- c(popcallList, list(site_color_category = ord_vars[1]))
	} 
	if( length(ord_vars) > 1){
		popcallList <- c(popcallList, list(site_shape_category = ord_vars[2]))
	}
	
	# Create the plot
	do.call("plot_ordination_phyloseq", popcallList)
}
################################################################################
################################################################################
#' Convert an otuTable object into a data.frame useful for plotting
#' in the ggplot2 framework.
#'
#' @usage otu2df(otu, taxavec, map, keepOnlyTheseTaxa=NULL, threshold=NULL)
#'
#' @param otu An \code{otuTable} object.
#'
#' @param taxavec A character vector of the desired taxonomic names to 
#'  categorize each species in otu. 
#' 
#' @param map The corresponding sampleData object for \code{otu}.
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
#' the sampleData component of \code{otu}.
#'
#' @usage taxaplot(otu, taxavec="Domain",
#'	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category,  
#'	facet_formula = . ~ TaxaGroup)
#'
#' @param otu (Required). An \code{otuTable} object, or higher-order object that contains
#'  an otuTable and sampleData (e.g. ``otuSam'' class and its superclasses.).
#'  If \code{otu} does not contain a taxTab slot (is a class that does not
#'  have ``Tax'' in its title), then the second argument, \code{taxavec}, is
#'  required and should have length equal to the number of species/taxa in
#'  \code{otu}.
#'
#' @param taxavec A character vector of the desired taxonomic names to 
#'  categorize each species in \code{otu}. If \code{otu} is a higher-order
#'  object that
#'  contains a taxonomyTable, then taxavec can alternatively specify the
#'  desired taxonomic level as a character string of length 1. 
#'  E.g. \code{taxavec = "Phylum"}. Default value is \code{"Domain"}.
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
#' @param x_category A character string indicating which sampleData column should be
#'  used to define the horizontal axis categories. Default is \code{"sample"}. Note
#'  that a few column-names are added by default and are available as options. 
#'  They are ``sample'', ``Abundance'', and ``TaxaGroup''.
#' 
#' @param fill_category A character string indicating which sampleData column
#'  should be used to define the fill color of the bars. This does not have to 
#'  match \code{x_category}, but does so by default. Note
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
#' @return A ggplot2 graphic object.
#'
#' @seealso \code{\link{otu2df}}, \code{\link{qplot}}, \code{\link{ggplot}}
#'
#' @export
#' @docType methods
#' @rdname taxaplot-methods
#'
#' @examples #
#' # data(ex1)
#' # taxaplot(ex1, "Class", threshold=0.85, x_category="Diet",
#' # fill_category="Diet", facet_formula = Gender ~ TaxaGroup)
setGeneric("taxaplot", function(otu, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category,  
	facet_formula = . ~ TaxaGroup){
		standardGeneric("taxaplot")
})
################################################################################
#' @rdname taxaplot-methods
#' @aliases taxaplot,phyloseq-method
setMethod("taxaplot", "phyloseq", function(otu, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category, 
	facet_formula = . ~ TaxaGroup){

	# Some preliminary assignments. Assumes otu has non-empty sampleData slot.
	map <- sampleData(otu)
	if( length(taxavec) == 1 ){ 
		taxavec <- as(taxTab(otu), "matrix")[, taxavec, drop=TRUE]
	}

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
			# data=df, eval(call("aes", x=as.name(x_category), 
				# y=quote(Abundance), fill=as.name(fill_category) )),
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
})
################################################################################
#' Create a taxa graph adjacency matrix from an \code{otuTable}.
#' 
#' By default, this uses a presence/absence criteria for taxa to determine
#' if samples (vertices) are connected by edges, and then plot the result
#' using the \code{\link[igraph]{igraph-package}}.
#'
#' @usage makenetwork(physeq, plotgraph=TRUE, 
#'			community=TRUE, threshold=0, incommon=0.4, method="jaccard")
#'
#' @param physeq (Required). An \code{\link{otuTable-class}}
#'  or \code{\link{phyloseq-class}} object.
#'
#' @param plotgraph A logical. Default TRUE. Indicates whether or 
#'  not to plot the network
#'
#' @param community A logical. Default TRUE. If TRUE, RETURN will
#'  include the community groups.
#'
#' @param threshold A single-value positive integer. Default 0. Indicates the
#'  number of individuals/observations required for presense to be counted.
#'  For example, threshold > 1 means that species with 2 or more reads
#'  are considered present.
#'
#' @param incommon The fraction of co-occurring samples required for a pair
#'  of species to be given an edge. Default is 0.4. This parameter has a 
#'  strong influence on the resulting network, and must be chosen carefully.
#'
#' @param method A character string indicating the default distance method
#'  parameter provided to \code{\link{vegdist}}. Default is ``jaccard''.
#'
#' @return Creates a plot of the adjacency graph and optionally returns the
#'  community groups.
#' 
#' @author Susan Holmes
#'
#' @export
#' @docType methods
#'
#' @examples
#' ## data(enterotype)
#' ## makenetwork(enterotype)
makenetwork <- function(physeq, plotgraph=TRUE, community=TRUE, threshold=0, incommon=0.4, method="jaccard"){	
	#require(vegan); require(igraph)

	abund <- otuTable(physeq)

	# transpose if speciesAreRows (vegan orientation)
	if( speciesAreRows(abund) ){ abund <- t(abund) }

	### Only take the rows where there are at least one value over threshold
	abundance <- abund[rowSums(abund) > threshold, ]
	n         <- nrow(abundance)

	# Convert to 1,0 binary matrix for input to vegdist. -0 converts to numeric
	presenceAbsence <- (abundance > threshold) - 0

	##Compute the Jaccard distance between the rows, this will only make points
	##closer if they are actually present together	
	##You could use any of the other distances in vegan or elsewhere
	jaccpa <- vegdist(presenceAbsence, method)
	###Distances in R are vectors by default, we make them into matrices	
	jaacm <- as.matrix(jaccpa)
	coinc <- matrix(0, n, n)
	ind1  <- which((jaacm>0 & jaacm<(1-incommon)),arr.ind=TRUE)
	coinc[ind1] <- 1
	dimnames(coinc) <- list(dimnames(abundance)[[1]],dimnames(abundance)[[1]])	
	###If using the network package create the graph with	
	###	g<-as.network.matrix(coinc,matrix.type="adjacency")
	####Here I use the igraph adjacency command
	ig=graph.adjacency(coinc)

	###Take out the isolates
	isolates=V(ig)[ degree(ig)==0 ]
	ignoisol=delete.vertices(ig, V(ig)[ degree(ig)==0 ])
	if (plotgraph==TRUE){
		plot(ignoisol, layout=layout.fruchterman.reingold, 
			vertex.size=0.6, vertex.label.dist=0.1, 
			edge.arrow.mode="-",vertex.color="red",
			vertex.label=NA,edge.color="blue")		
		title("Co-occurrence graph without isolates")
	}
    if (community==TRUE){
		communitywalk=walktrap.community(ignoisol)
		nonisolates= V(ig)[ degree(ig)!=0 ]
		group0=nonisolates[which(communitywalk$membership==0)]
		group1=nonisolates[which(communitywalk$membership==1)]
		###You can then play around with coloring and labelling in the graph
		###For help don't forget to look up plot.igraph or plot.network
		###not just plot as it inherits the plot method appropriate to its class
		groups=list(group0,group1)
		return(groups)
	}
}
################################################################################
# tipsymbols and tiptext. Need to make both tipsymbols and tiptext documentation
# point to this one.
################################################################################
#' Annotate tips on a tree with symbols or text.
#'
#' There were some unexpected behavior from the \code{\link[ape]{tiplabels}}
#' function in ape. These
#' functions are intended to act as simplified versions that act as a convenience
#' wrapper for \code{points()} or \code{text()} functions, respectively, 
#' but where the tip
#' coordinates are specified by giving the tip ID (integer) as input.
#' For \code{tiptext()}, make sure to include a \code{labels=} argument, which
#' will be passed on to \code{\link[graphics]{text}}.
#'
#' @usage tipsymbols(tip, adj=c(0.5, 0.5), ...)
#' @usage tiptext(tip, adj=c(0.5, 0.5), ...)
#'
#' @param tip An integer specifying the tip ID in a tree that for which the 
#'  base plot has already been generated and is still available to \code{R}.
#' 
#' @param adj A 2 element numeric vector specifying a position adjustment.
#' 
#' @param ... Additional plotting parameters that are passed to 
#'  \code{\link{points}} or \code{\link{text}} in the R base graphics.
#'  Again, for \code{tiptext()}, make sure to include a \code{labels=} argument.
#'
#' @return No objects returned. Symbol or text is plotted on the available
#'  graphic device.
#'
#' @export
#' @rdname tip-annotate
#' 
#' @seealso \code{\link[ape]{tiplabels}}, \code{\link[graphics]{points}}, \code{\link[graphics]{text}}
#' @examples #
#' ## data(ex1)
#' ## # for reproducibility
#' ## set.seed(711)
#' ## ex2 <- prune_species(sample(species.names(ex1), 50), ex1)
#' ## plot( tre(ex2) )
#' ## tipsymbols(pch=19)
#' ## tipsymbols(1, pch=22, cex=3, col="red", bg="blue")
#' ## tiptext(2, labels="my.label")
tipsymbols <- function(tip, adj=c(0.5, 0.5), ...){
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if ( missing(tip) ){ 
		tip <- 1:lastPP$Ntip
	}
    XX <- lastPP$xx[tip]
    YY <- lastPP$yy[tip]
	points( (XX + adj[1] - 0.5), (YY + adj[2] - 0.5), ... )
}
################################################################################
# Custom text plotting function
#' @export
#' @rdname tip-annotate
tiptext <- function(tip, adj=c(0.5, 0.5), ...){
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if ( missing(tip) ){ 
		tip <- 1:lastPP$Ntip
	}
    XX <- lastPP$xx[tip]
    YY <- lastPP$yy[tip]
	text( (XX + adj[1] - 0.5), (YY + adj[2] - 0.5), ... )
}
################################################################################
#' Plot tree with easy tip annotation.
#'
#' Requires a \code{\link{phyloseq-class}} that contains a tree (\code{\link{tre}}), 
#' sample data (\code{\link{sampleData}}),
#' and abundance table (\code{\link{otuTable}}).
#'
#' @usage plot_tree_phyloseq(physeq, color_factor=NULL, shape_factor=NULL, 
#'  base_size=1, size_scaling_factor = 0.2, opacity=2/3,
#'  custom_color_scale=NULL, custom_shape_scale=NULL, 
#'  type_abundance_value=FALSE, printTheseTaxa=NULL, treeTitle="Annotated Tree", ...)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} with non-empty 
#'  tree, sampleData, and otuTable components.
#'
#' @param color_factor A character string specifying the column
#'  of the sampleData that will be used for setting the color of symbols.
#'
#' @param shape_factor A character string specifying the column
#'  of the sampleData that will be used for setting the shape of symbols.
#'  
#' @param base_size The minimum size expansion factor of symbols plotted next
#'  to tips. The default value is 1. 
#'
#' @param size_scaling_factor A numeric, greater than or equal to 0, that is
#'  multiplied by the log10 of taxa abundance; the product of which is summed
#'  with the \code{base_size} argument to determine the size scaling factor provided
#'  to \code{\link{tipsymbols}}. The default value is 0.15. The larger the value,
#'  the larger the symbols representing sites with many idividuals of a
#'  particular taxa. A value of zero means there will be no scaling symbol
#'  size by the abundance value.
#'
#' @param opacity The opacity (or alpha value). Numeric between 0, 1.
#'  Defaul value is 2/3.
#'
#' @param custom_color_scale A character vector of the desired custom color scale.
#'  This should
#'  be a scale, not an aesthetic map. Therefore, it will in most-cases 
#'  contain only unique elements, unless two different categories of data
#'  are supposed to have the same color. Default value is NULL, which
#'  invokes a default color scale using the \code{\link{rainbow}} function.
#' 
#' @param custom_shape_scale An integer vector of values in the categorical
#'  scale of symbol shapes, analogous to \code{custom_color_scale}. Default
#'  value is \code{NULL}, which uses the fill-able symbols described in
#'  \code{\link{points}}, beginning with 21. 
#' 
#' @param type_abundance_value Logical. Whether or not the otuTable value
#'  (the number of individuals, typically) should be added to the center
#'  of symbols when the value is greater than one. 
#'  Default is FALSE, indicating no labels.
#'
#' @param printTheseTaxa a character vector of the taxa names in \code{physeq}
#'  that should be labeled on the tree plot adjacent to the right. Default is
#'  NULL. Not yet implemented.
#'
#' @param treeTitle (Optional). Character string, for the title
#'  of the graphic. Default is \code{"Annotated Tree"}.
#'
#' @param ... Additional parameters passed on to \code{\link{tipsymbols}}.
#'
#' @return Creates a phylogenetic tree, with additional symbols annotated on
#'  each tip to indicate in which samples the particular taxa was observed.
#'
#' @export
#'
#' @examples
#' # data(ex1)
#' # ex2 <- ex1
#' ## # for reproducibility
#' ## set.seed(711)
#' # species.names(ex2) <- sample(species.names(ex1), 50)
#' # plot_tree_phyloseq(ex2)
#' # plot_tree_phyloseq(ex2, shape_factor="Diet")
#' # plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet")
#' # plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet", 
#' 	# size_scaling_factor=0.6, type_abundance_value=TRUE)
#' # plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet", 
#'	# size_scaling_factor=0.6, custom_color_scale=c("blue", "magenta"))
#' # plot(phyloseqTree(ex2), color_factor="Gender", shape_factor="Diet", 
#'  # size_scaling_factor=0.6, custom_color_scale=c("blue", "magenta") )
plot_tree_phyloseq <- function(physeq, color_factor=NULL, shape_factor=NULL, 
	base_size=1, size_scaling_factor = 0.2, opacity=2/3,
	custom_color_scale=NULL, custom_shape_scale=NULL, 
	type_abundance_value=FALSE, printTheseTaxa=NULL, treeTitle="Annotated Tree", ...){

	# Initialize categories vector, giving the sampleData index of aesthetics.
	categories <- character(2); names(categories) <- c("color", "shape")
	
	# Only look for default factors if color_factor AND shape_factor absent.
	if( is.null(color_factor) & is.null(shape_factor) ){
		# Want to use up to the first two factors in sampleData as default.
		smdf_factor_types <- lapply(data.frame(sampleData(physeq)), class) == "factor"
		available_factors <- colnames(sampleData(physeq))[	smdf_factor_types ]
		categories["color"] <- available_factors[1]
		if( length(available_factors) > 1 ){
			categories["shape"] <- available_factors[2]
		}	
	}
	
	# If categories were specified by user, get them.
	if(!is.null(color_factor)){ categories["color"] <- color_factor }
	if(!is.null(shape_factor)){ categories["shape"] <- shape_factor }
	
	# color
	if( categories["color"] != "" ){
		color_factor <- data.frame(sampleData(physeq))[, categories["color"]]
		# determine the color scale.
		if( is.null(custom_color_scale) ){
			avail_colors <- rainbow(length(levels(color_factor)), alpha=opacity)
		} else {
			avail_colors <- custom_color_scale
		}
		names(avail_colors) <- levels(color_factor) 
		# set the color vector.
		color_vec    <- avail_colors[color_factor]
	} else {
		color_vec <- rep("black", nsamples(physeq))
	}
	names(color_vec) <- sample.names(physeq)

	# shape	 custom_shape_scale
	if( categories["shape"] != "" ){
		shape_factor <- data.frame(sampleData(physeq))[, categories["shape"]]
		# determine the shape scale.
		if( is.null(custom_shape_scale) ){
			avail_shapes <- (1:(length(levels(shape_factor)))) + 20
		} else {
			avail_shapes <- custom_shape_scale
		}
		names(avail_shapes) <- levels(shape_factor) 
		# set the shape vector
		shape_vec <- avail_shapes[shape_factor]
	} else {
		shape_vec <- rep(22, nsamples(physeq))
	}
	names(shape_vec) <- sample.names(physeq)
	
	# Access phylo-class tree. Error if it is missing.
	tree <- tre(physeq)
	
	# Now plot the initial, symbol-less tree. Must be first to get the proper
	# x, y limits to calculate the scales of the annotation objects.
	ape::plot.phylo(tree, type="phylogram", show.tip.label=FALSE,
		xpd=NA, no.margin=TRUE, edge.width=1.5)

	# Store information about the tree-plot
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	#str(lastPP)
	xlims <- lastPP$x.lim
	ylims <- lastPP$y.lim	

	# Add title of tree
	text(x=(0.35*xlims[2]), y=(1.02*ylims[2]), labels=treeTitle[1], pos=4, cex=1.5)

	# Add scale bar
	add.scale.bar(x=(0.5*xlims[2]), y=0, length=0.1)
	adj.j.start <- 0.5 + (0.01 * xlims[2])
	adj.step    <- 0.0155 * xlims[2]
	########################################
	# Now loop over sample variables
	########################################
	for( i in species.names(physeq) ){ # i loops over species
		#i = species.names(physeq)[speciesSums(physeq) == max(speciesSums(physeq))]
		# sitesi should hold the sites that are relevant to species "i"
		sitesi <- getSamples( physeq, i)
		sitesi <- sitesi[sitesi >= 1] 
		
		# assign to pchi / coli values if species "i" is in the associated site
		tipi  <- which( tree$tip.label %in% i )
		
		# Now loop through to indicate multiple symbols if multiple samples
		adj.j <- adj.j.start
		for( j in names(sitesi) ){ # j loops over sites w/in a given species
			size_i = base_size + size_scaling_factor * log(sitesi[j], 10)
			tipsymbols(tip=tipi, adj=c(adj.j, 0.5), pch=shape_vec[j], bg=color_vec[j],
				col=color_vec[j], cex=size_i, ...)
				
			# print the number of individuals observed of this species if > 1
			# And if user has asked to plot the abundance values.
			if( type_abundance_value ){
				if( sitesi[j] > 1 ){
					tiptext(tip=tipi, adj=c(adj.j, 0.5),
						labels=as.character(sum(sitesi[j])),
						col="black", cex=(size_i/3))
				}				
			}
			
			# Increase the horizontal coordinate by adj.step, for next loop.
			adj.j <- adj.j + adj.step
		}
		
		# Add taxa names if requested.
		# if( i %in% printTheseTaxa ){
			# ref_species_name <- species_table[i, "taxa"]
			# tiptext(tip=tipi, adj=c(adj.j-0.03, 0.3), labels=ref_species_name,
				# col="black", pos=4, cex=(cex.symbol/3))		
		# }	
	}
	
	# Need to add legends module. Probably best accomplished with a separate margin
	# on the right side of the plot.	

}
################################################################################
# colname2hex
################################################################################
#' Convert color name to hex value.
#'
#' For converting color name to color hex and providing an optional
#' opacity parameter (alpha)
#'
#' @param colname A character vector of \code{R} color names.
#' 
#' @param alpha The standard alpha parameter specifying opacity/transparency.
#'
#' @return A color hex value of length equal to length of \code{colname}.
#' @seealso col2rgb rgb2hsv colors
#'
#' @keywords internal 
colname2hex <- function(colname, alpha=1){
	if( length(colname) == 1){
		do.call("hsv", c(as.list(rgb2hsv(col2rgb(colname))[, 1]), alpha=alpha))
	} else if( length(colname) >1 ){
		sapply(colname, colname2hex, alpha)
	}	
}
################################################################################