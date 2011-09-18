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
#' @param mod A \code{cca} or \code{rda} results object. 
#'  See \code{\link{cca.object}}. For phyloseq
#'  objects, you probably want to see \code{\link{cca.phyloseq}} or 
#'  \code{\link{rda.phyloseq}}.
#'  
#' @param object A \code{phyloseq} object or one of its child classes
#'  (superclasses). Required. Necessary because ordination result objects
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
#'  about the samples/sites that is not already available in the \code{sampleMap} of
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
#' @return A ggplot2 plot object.
#'
#' @seealso plot.cca cca.object
#' @import vegan
#' @export
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
	# if there is a sampleMap, cbind it to the sites data
	if( !is.null(sampleMap(object)) ){
		sitesdata <- cbind(sitesdata, data.frame(sampleMap(object)))
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
#' @param X A formula object, with the left-hand side specifying a single 
#'  phyloseq object that contains (at minimum) an \code{otuTable} and a 
#'  \code{sampleMap}. The right-hand side expects at least two different
#'  variates, present in the \code{sampleMap} of the LHS. For available
#'  variate names in your object, try \code{colnames(sampleMap(ex1))}, where
#'  \code{ex1} is the phyloseq object containing your data. Because this is
#'  a formula object, quotes should not be used. See \code{\link{formula}}
#'  for details about writing a formula in R.
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
calcplot <- function(X, RDA_or_CCA="cca", object=get(all.vars(X)[1]), ...){
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

	# default is to use the R.H.S. elements of X as color and shape
	plot_ordination_phyloseq(mod,
		object = object, 
		site_shape_category = ord_vars[1],
		site_color_category = ord_vars[2]
	)
}
################################################################################
