#
# extension of plot methods for phyloseq object.
# 
################################################################################
################################################################################
################################################################################
################################################################################
#' Generic plot defaults for phyloseq.
#'
#' There are many useful examples of phyloseq graphics functions in the
#' \href{http://joey711.github.com/phyloseq/}{phyloseq online tutorials}.
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
#'  \href{https://github.com/joey711/phyloseq/wiki/Graphics-Examples}{phyloseq graphics examples (wiki)}.
#'
#'  \code{\link{plot_ordination}}
#'  \code{\link{plot_heatmap}}
#'  \code{\link{plot_tree}}
#'  \code{\link{plot_network}}
#'  \code{\link{plot_taxa_bar}}
#'  \code{\link{plot_richness}}
#'
#' @export
#' @docType methods
#' @rdname plot_phyloseq-methods
#'
#' @examples 
#' data(esophagus)
#' plot_phyloseq(esophagus)
setGeneric("plot_phyloseq", function(physeq, ...){ standardGeneric("plot_phyloseq") })
#' @aliases plot_phyloseq,phyloseq-method
#' @rdname plot_phyloseq-methods
setMethod("plot_phyloseq", "phyloseq", function(physeq, ...){
	if( all(c("otu_table", "sample_data", "phy_tree") %in% getslots.phyloseq(physeq)) ){
		plot_tree(esophagus, color="samples")	
	} else if( all(c("otu_table", "sample_data", "tax_table") %in% getslots.phyloseq(physeq) ) ){
		plot_bar(physeq, ...)
	} else if( all(c("otu_table", "phy_tree") %in% getslots.phyloseq(physeq)) ){
		plot_tree(esophagus, color="samples")	
	} else {
		plot_richness(physeq)
	}
})
################################################################################
################################################################################
#' Plot a network using ggplot2 (represent microbiome)
#'
#' There are many useful examples of phyloseq network graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_network-examples}{phyloseq online tutorials}.
#' A custom plotting function for displaying networks
#' using advanced \code{\link[ggplot2]{ggplot}}2 formatting.
#' The network itself should be represented using
#' the \code{igraph0} package.
#' For the \code{\link{phyloseq-package}} it is suggested that the network object
#' (argument \code{g})
#' be created using the
#'  \code{\link{make_network}} function, 
#' and based upon sample-wise or taxa-wise microbiome ecological distances 
#' calculated from a phylogenetic sequencing experiment 
#' (\code{\link{phyloseq-class}}).
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage plot_network(g, physeq=NULL, type="samples", 
#' 	color=NULL, shape=NULL, point_size=4, alpha=1,
#' 	label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
#'
#' @param g (Required). An \code{igraph0}-class object created
#'  either by the convenience wrapper \code{\link{make_network}}, 
#'  or directly by the tools in the igraph0-package.
#'
#' @param physeq (Optional). Default \code{NULL}. 
#'  A \code{\link{phyloseq-class}} object on which \code{g} is based.
#'
#' @param type (Optional). Default \code{"samples"}.
#'  Whether the network represented in the primary argument, \code{g},
#'  is samples or taxa/OTUs.
#'  Supported arguments are \code{"samples"}, \code{"taxa"},
#'  where \code{"taxa"} indicates using the taxa indices,
#'  whether they actually represent species or some other taxonomic rank.
#'
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#' 
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).
#' 
#' @param point_size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#' 
#' @param alpha (Optional). Default \code{1}.
#'  A value between 0 and 1 for the alpha transparency of the vertex points.
#' 
#' @param label (Optional). Default \code{"value"}.
#'  The name of the sample variable in \code{physeq} to use for 
#'  labelling the vertex points.
#' 
#' @param hjust (Optional). Default \code{1.35}.
#'  The amount of horizontal justification to use for each label.
#' 
#' @param line_weight (Optional). Default \code{0.3}.
#'  The line thickness to use to label graph edges.
#' 
#' @param line_color (Optional). Default \code{color}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of lines (graph edges).
#' 
#' @param line_alpha (Optional). Default \code{0.4}.
#'  The transparency level for graph-edge lines.
#'
#' @param layout.method (Optional). Default \code{layout.fruchterman.reingold}.
#'  A function (closure) that determines the placement of the vertices
#'  for drawing a graph. Should be able to take an \code{igraph0}-class
#'  as sole argument, and return a two-column coordinate matrix with \code{nrow}
#'  equal to the number of vertices. For possible options already included in 
#'  \code{igraph0}-package, see the others also described in the help file:
#' 
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' \code{\link[igraph0]{layout.fruchterman.reingold}}
#'
#' @return A \code{\link{ggplot}}2 plot representing the network,
#'  with optional mapping of variable(s) to point color or shape.
#' 
#' @seealso 
#'  \code{\link{make_network}}
#'
#' @references
#'  This code was adapted from a repo original hosted on GitHub by Scott Chamberlain:
#'  \url{https://github.com/SChamberlain/gggraph}
#'
#'  The code most directly used/modified was first posted here:
#'  \url{http://www.r-bloggers.com/basic-ggplot2-network-graphs/}
#' 
#' @import ggplot2
#' @import reshape
#' @importFrom igraph0 layout.fruchterman.reingold
#' @importFrom igraph0 get.edgelist
#' @export
#' @examples 
#' 
#' data(enterotype)
#' ig <- make_network(enterotype, max.dist=0.3)
#' plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # Change distance parameter
#' ig <- make_network(enterotype, max.dist=0.2)
#' plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
plot_network <- function(g, physeq=NULL, type="samples", 
	color=NULL, shape=NULL, point_size=4, alpha=1,
	label="value", hjust = 1.35, 
	line_weight=0.5, line_color=color, line_alpha=0.4,
	layout.method=layout.fruchterman.reingold, title=NULL){

	# disambiguate species/OTU/taxa as argument type...
	if( type %in% c("taxa", "species", "OTUs", "otus", "otu") ){
		type <- "taxa"
	}

	# Make the edge-coordinates data.frame
	edgeDF    <- data.frame(get.edgelist(g))
	edgeDF$id <- 1:length(edgeDF[, 1])

	# Make the vertices-coordinates data.frame
	vertDF    <- layout.method(g)
	colnames(vertDF) <- c("x", "y")
	vertDF    <- data.frame(value=g[[9]][[3]][["name"]], vertDF)
	
	# If phyloseq object provided,
	# AND it has the relevant additional data
	# THEN add it to vertDF
	if( !is.null(physeq) ){
		extraData <- NULL
		if( type == "samples" & !is.null(sample_data(physeq, FALSE)) ){
			extraData <- sample_data(physeq)[as.character(vertDF$value), ]
		} else if( type == "taxa" & !is.null(tax_table(physeq, FALSE)) ){
			extraData <- tax_table(physeq)[as.character(vertDF$value), ]
		}
		# Only mod vertDF if extraData exists
		if( !is.null(extraData) ){
			vertDF <- data.frame(vertDF, extraData) 			
		}
	}

	# Combine vertex and edge coordinate data.frames
	graphDF   <- merge(melt(edgeDF, id="id"), vertDF, by = "value") 
 
	# Initialize the ggplot
	p <- ggplot(vertDF, aes(x, y)) 

	# Strip all the typical annotations from the plot, leave the legend
	p <- p + theme_bw() + 
			theme(
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), 
				axis.text.x      = element_blank(),
				axis.text.y      = element_blank(),
				axis.title.x     = element_blank(),
				axis.title.y     = element_blank(),
				axis.ticks       = element_blank(),
				panel.border     = element_blank()
			)

	# Add the graph vertices as points
	p <- p + geom_point(aes_string(color=color, shape=shape), size=point_size, na.rm=TRUE)

	# Add the text labels
	if( !is.null(label) ){
		p <- p + geom_text(aes_string(label=label), size = 2, hjust=hjust, na.rm=TRUE)
	}
	
	# Add the edges:
	p <- p + geom_line(aes_string(group="id", color=line_color), 
				graphDF, size=line_weight, alpha=line_alpha, na.rm=TRUE)
				
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	return(p)
}
################################################################################
################################################################################
#' Plot richness estimates, flexibly with ggplot2
#'
#' There are many useful examples of phyloseq richness graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_richness-examples}{phyloseq online tutorials}.
#' Performs a number of standard richness estimates using the 
#' \code{\link{estimate_richness}} function,
#' and returns a \code{ggplot} plotting object. 
#' This plot shows the individual richness estimates for each
#' sample, as well as the observed richness. 
#' You must use untrimmed datasets
#' for meaningful results, as these estimates (and even the ``observed'' richness)
#' are highly dependent on the number of singletons. You can always trim the data
#' later on if needed, just not before using this function.
#'
#'  NOTE: Because this plotting function incorporates the output from 
#'  \code{\link{estimate_richness}}, the variable names of that output should
#'  not be used as \code{x} or \code{color} (even if it works, the resulting
#'  plot might be kindof strange, and not the intended behavior of this function).
#'  The following are the names you will want to avoid using in \code{x} or \code{color}:
#'
#'  \code{c("S.obs", "S.chao1", "se.chao1", "S.ACE", "se.ACE", "shannon", "simpson")}
#'
#' @usage plot_richness(physeq, x="samples", color=NULL, shape=NULL, title=NULL, shsi=FALSE)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otu_table-class}}. The data about which you want to estimate
#'  the richness.
#'
#' @param x (Optional). A variable to map to the horizontal axis. The vertical
#'  axis will be mapped to richness estimates and have units of total taxa.
#'  This parameter (\code{x}) can be either a character string indicating a
#'  variable in \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'
#'  The default value is \code{"samples"}, which will map each sample's name
#'  to a separate horizontal position in the plot.
#'
#' @param color (Optional). Default \code{NULL}. The sample variable to map
#'  to different colors. Like \code{x}, this can be a single character string 
#'  of the variable name in 
#'  \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'  The color scheme is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_color_manual}}.
#'
#' @param shape (Optional). Default \code{NULL}. The sample variable to map
#'  to different shapes. Like \code{x} and \code{color},
#'  this can be a single character string 
#'  of the variable name in 
#'  \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'  The shape scale is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_shape_manual}}.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @param shsi (Optional). Default \code{FALSE}. Logical.
#'  Whether or not to include Shannon and Simpson indices
#'  in the graphic as well.
#'
#' @return A \code{\link{ggplot}} plot object summarizing
#'  the richness estimates, and their standard error.
#' 
#' @seealso 
#'  \code{\link{estimate_richness}},
#'  \code{\link[vegan]{estimateR}},
#'  \code{\link[vegan]{diversity}}
#'
#' There are many more interesting examples at the
#' \href{http://joey711.github.com/phyloseq/plot_richness-examples}{phyloseq online tutorials}.
#'
#' @import ggplot2
#' @import reshape
#' @export
#' @examples 
#' ## There are many more interesting examples at the phyloseq online tutorials.
#' ## http://joey711.github.com/phyloseq/plot_richness-examples
#' data(GlobalPatterns)
#' GP = prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
#' plot_richness(GP, x = "SampleType", color="SampleType")
#' plot_richness(GP, x = "SampleType", color="SampleType", shsi=TRUE)
plot_richness <- function(physeq, x="samples", color=NULL, shape=NULL, title=NULL, shsi=FALSE){

	# Make the plotting data.frame 
	DF <- data.frame(estimate_richness(physeq), sample_data(physeq))
	
	# If there is no "samples" variable in DF, add it
	if( !"samples" %in% names(DF) ){
		DF$samples = sample_names(physeq)
	}
	
	# sample_names used to be default, and should also work.
	# #backwardcompatibility
	if( !is.null(x) ){
		if( x %in% c("sample", "samples", "sample_names") ){
			x = "samples"
		}
	}

	# Define "measure" variables and s.e. labels (ses).
	measures = c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson")
	ses = c("se.obs", "se.chao1", "se.ACE", "se.shannon", "se.simpson")

	# melt, for different richnesses...
	mdf = melt(DF, measure.vars=measures)

	# Merge s.e. into one "se" column
	mdf$se = NA_integer_
	mdf$wse = paste("se.", substr(mdf$variable, 3, 10), sep="")
	for( i in 1:nrow(mdf) ){
		if( mdf[i, "wse"] %in% c("se.chao1", "se.ACE") ){
			mdf[i, "se"] = mdf[i, (mdf[i, "wse"])]
		}
	}
	
	# Rm shannon/simpson if !shsi
	if( !shsi ){
		mdf = subset(mdf, variable %in% measures[1:3])
	}
	
	# map variables
	richness_map <- aes_string(x=x, y="value", color=color, shape=shape)		
	
	# Make the ggplot.
	p <- ggplot(mdf, richness_map) + 
		geom_point(na.rm=TRUE) + 
		geom_errorbar(aes(ymax=value + se, ymin=value - se), width=0.2) +	
		theme(axis.text.x = element_text(angle = -90, hjust = 0))
	
	# Add label according to whether or not shannon/simpson indices are included
	if(shsi){
		p = p + ylab('Alpha Diversity Measure') 				
	} else {
		p = p + ylab('Richness [number of taxa]') 		
	}
		
	# Facet differently, depending on whether shannon or simpson indices included
	if(shsi){
		p = p + facet_wrap(~variable, nrow=1, scales="free") 				
	} else {
		p = p + facet_grid(~variable) 		
	}
		
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	return(p)
}
################################################################################
################################################################################
# The general case, could plot samples, taxa, or both (biplot/split). Default samples.
################################################################################
#' General ordination plotter based on ggplot2.
#'
#' There are many useful examples of phyloseq ordination graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_ordination-examples}{phyloseq online tutorials}.
#' Convenience wrapper for plotting ordination results as a 
#' \code{ggplot2}-graphic, including
#' additional annotation in the form of shading, shape, and/or labels of
#' sample variables.
#'
#' @usage plot_ordination(physeq, ordination, type="samples", axes=c(1, 2),
#'	color=NULL, shape=NULL, label=NULL, title=NULL, justDF=FALSE)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}. 
#'  The data about which you want to 
#'  plot and annotate the ordination.
#'
#' @param ordination (Required). An ordination object. Many different classes
#'  of ordination are defined by \code{R} packages. Ordination classes
#'  currently supported/created by the \code{\link{ordinate}} function are
#'  supported here. There is no default, as the expectation is that the 
#'  ordination will be performed and saved prior to calling this plot function.
#'
#' @param type (Optional). The plot type. Default is \code{"samples"}. The
#'  currently supported options are 
#'  \code{c("samples", "sites", "species", "taxa", "biplot", "split", "scree")}.
#'  The option
#'  ``taxa'' is equivalent to ``species'' in this case, and similarly,
#'  ``samples'' is equivalent to ``sites''. 
#'  The options
#'  \code{"sites"} and \code{"species"} result in a single-plot of just the 
#'  sites/samples or species/taxa of the ordination, respectively.
#'  The \code{"biplot"} and \code{"split"} options result in a combined
#'  plot with both taxa and samples, either combined into one plot (``biplot'')
#'  or 
#'  separated in two facet panels (``split''), respectively.
#'  The \code{"scree"} option results in a call to \code{\link{plot_scree}},
#'  which produces an ordered bar plot of the normalized eigenvalues
#'  associated with each ordination axis. 
#'
#' @param axes (Optional). A 2-element vector indicating the axes of the 
#'  ordination that should be used for plotting. 
#'  Can be \code{\link{character-class}} or \code{\link{integer-class}},
#'  naming the index name or index of the desired axis for the horizontal 
#'  and vertical axes, respectively, in that order. The default value, 
#'  \code{c(1, 2)}, specifies the first two axes of the provided ordination.
#'
#' @param color (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map to
#'  colors in the plot. 
#'  This can be a sample variable 
#'  (among the set returned by \code{sample_variables(physeq)} )
#'  or
#'  taxonomic rank
#'  (among the set returned by \code{rank_names(physeq)}).
#'  
#'  Alternatively, if \code{type} indicates a single-plot 
#'  (\code{"samples"} or \code{"species"}), then
#'  it is also possible to supply a custom vector with length equal to
#'  the relevant number of samples or species
#'  (\code{nsamples(physeq)} or \code{ntaxa(physeq)}).
#' 
#'  Finally,
#'  The color scheme is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_color_manual}}.
#'
#' @param shape (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map
#'  to different shapes on the plot. 
#'  Similar to \code{color} option, but for the shape if points.
#' 
#'  The shape scale is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_shape_manual}}.
#'
#' @param label (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map to text labels on the plot.
#'  Similar to \code{color} option, but for plotting text.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @param justDF (Optional). Default \code{FALSE}. Logical.
#'  Instead of returning a ggplot2-object, do you just want the relevant
#'  \code{data.frame} that was used to build the plot? This is a 
#'  user-accessible option for obtaining the \code{data.frame}, in 
#'  in principal to make a custom plot that isn't possible with the
#'  available options in this function. For contributing new functions
#'  (developers), the  
#'  \code{\link{phyloseq-package}} provides/uses an internal function
#'  to build the key features of the \code{data.frame} prior to plot-build.
#'
#' @return A \code{\link{ggplot}} plot object, graphically summarizing
#'  the ordination result for the specified axes.
#' 
#' @seealso 
#'  The examples on the phyloseq wiki page for \code{plot_ordination} show 
#'  many more examples:
#'
#' \url{https://github.com/joey711/phyloseq/wiki/plot_ordination}
#'
#' Also see the general wrapping function:
#'
#' \code{\link{plot_phyloseq}}
#'
#' @import ggplot2
#' @export
#' @examples 
#' data(GlobalPatterns)
#' # Need to clean the zeros from GlobalPatterns:
#' GP <- prune_taxa(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
#' # Define a human-associated versus non-human binary variable:
#' sample_data(GP)$human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
#' # Get the names of the most-abundant
#' top.TaxaGroup <- sort(
#'    tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum, na.rm = TRUE),
#'    decreasing = TRUE)
#' top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
#' # Now prune further, to just the most-abundant phyla
#' GP <- subset_taxa(GP, Phylum %in% names(top.TaxaGroup))
#' topsp <- names(sort(taxa_sums(GP), TRUE)[1:200])
#' GP1   <- prune_taxa(topsp, GP)
#' GP.dpcoa <- ordinate(GP1, "DPCoA")
#' plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum")
#' # Customize with ggplot2 layers added directly to output
#' library("ggplot2")
#' plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType") + geom_line() + geom_point(size=5)
#' p <- plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType", shape="human")
#' print(p)
#' # library("ggplot2")
#' # p + geom_line() + geom_point(size=5)
#' # plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum") + geom_line() + geom_point(size=5)
#' plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum", label="SampleType")
#' plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum")
#' plot_ordination(GP1, GP.dpcoa, type="biplot", color="Phylum")
#' plot_ordination(GP1, GP.dpcoa, type="biplot", label="Phylum")
#' plot_ordination(GP1, GP.dpcoa, type="split", color="Phylum", label="SampleType")
#' plot_ordination(GP1, GP.dpcoa, type="split", color="SampleType", shape="Phylum", label="SampleType")
#' plot_ordination(GP1, GP.dpcoa, type="scree")
plot_ordination <- function(physeq, ordination, type="samples", axes=c(1, 2),
	color=NULL, shape=NULL, label=NULL, title=NULL, justDF=FALSE){

	if(class(physeq)!="phyloseq"){
		warning("Full functionality requires physeq be phyloseq-class with multiple components.")
	}
	official_types = c("sites", "species", "biplot", "split", "scree")
	if(type == "samples"){type <- "sites"} # vegan compatibility with phyloseq
	if(type == "taxa"){type <- "species"} # vegan compatibility with phyloseq
	if( !type %in% official_types ){
		warning("type argument not supported. type set to \"samples\".")
		type = "sites"
	}
	# Stop early by passing to plot_scree() if "scree" was chosen as a type
	if( type %in% c("scree") ){
		return( plot_scree(ordination, title=title) )
	}

	# Build data.frame:
	if( type %in% c("sites", "species") ){
		# Because of the way scores()/coord are bound first in DF, the first two axes should
		# always be x and y, respectively. This is also contingent on the "choices" argument
		# to scores() working properly		
		DF <- ord.plot.DF.internal(physeq, ordination, type, axes)
		# Add, any custom-supplied plot-mapped variables
		if( length(color) > 1 ){
			DF$color <- color
			names(DF)[names(DF)=="color"] <- deparse(substitute(color))
			color <- deparse(substitute(color))
		}
		if( length(shape) > 1 ){
			DF$shape <- shape
			names(DF)[names(DF)=="shape"] <- deparse(substitute(shape))
			shape <- deparse(substitute(shape))
		}	
		if( length(label) > 1 ){
			DF$label <- label
			names(DF)[names(DF)=="label"] <- deparse(substitute(label))
			label <- deparse(substitute(label))
		}
		x <- names(DF)[1]
		y <- names(DF)[2]			
	} else if( type %in% c("split", "biplot") ){
		# Define DFs
		specDF <- ord.plot.DF.internal(physeq, ordination, type="species", axes)
		siteDF <- ord.plot.DF.internal(physeq, ordination, type="sites", axes)
		# Define x-label and y-label before merge, use sample-axis names (arbitrary)
		names(specDF)[1] <- x <- names(siteDF)[1] # "x-axis"
		names(specDF)[2] <- y <- names(siteDF)[2] # "y-axis"
		# Add id.type label
		specDF$id.type <- "taxa"
		siteDF$id.type <- "samples"
		# Merge the two data.frame together, for joint plotting.
		DF <- merge(specDF, siteDF, all=TRUE)
		# Replace NA with "samples" or "taxa", where appropriate (factor/character)
		if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "samples") }
		if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "taxa") }
		if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "samples") }
		if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "taxa") }		
	}
	
	# In case user wants the plot-DF for some other purpose, return early
	if(justDF){return(DF)}

	# If there is nothing to map (data.frame only has two columns), just return simple plot
	if(ncol(DF)<=2){
		ord_map <- aes_string(x=x, y=y)
		p <- ggplot(DF, ord_map) + geom_point(na.rm=TRUE)
		return(p)
	}
	
	# Mapping section
	if( type %in% c("sites", "species", "split") ){
		ord_map <- aes_string(x=x, y=y, color=color, shape=shape, na.rm=TRUE)
	} else if(type=="biplot"){
		# biplot, id.type must map to color or size. Only color if none specified.
		if( is.null(color) ){
			ord_map <- aes_string(x=x, y=y, color="id.type",
							shape=shape, na.rm=TRUE)
		} else {
			ord_map <- aes_string(x=x, y=y, size="id.type",
							color=color, shape=shape, na.rm=TRUE)
		}
	}

	# Plot-building section
	p <- ggplot(DF, ord_map) + geom_point(na.rm=TRUE)
	
	# split/facet color and shape can be anything in one or other.
	if( type=="split" ){
		# split-option requires a facet_wrap
		p <- p + facet_wrap(~id.type, nrow=1)
	}
	
	# If biplot, adjust scales
	if( type=="biplot" ){	
		if( is.null(color) ){
			# Rename color title in legend.
			p <- update_labels(p, list(colour = "type")) #p + scale_color_discrete(name="type")
		} else {
			# Check if variable is discrete. 
			if( is.discrete(DF[, color]) ){
				# The following function reproduces ggplot2's default color scale.
				# From: http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
				gg_color_hue <- function(n) {
					hues = seq(15, 375, length=n+1)
					hcl(h=hues, l=65, c=100)[1:n]
				}
				colvals <- gg_color_hue(length(levels(as(DF[, color], "factor"))))
				names(colvals) <- levels(as(DF[, color], "factor"))
				# Now make the taxa or samples dark grey
				colvals[names(colvals) %in% c("samples", "taxa")] <- "grey45"
				# Now add the manually re-scaled layer with taxa/samples as grey
				p <- p + scale_colour_manual(values=colvals)
			}
			# Adjust size so that samples are bigger than taxa by default.
			p <- p + scale_size_manual("type", values=c(samples=5, taxa=2))		
		}
	}

	# Add the text labels
	if( !is.null(label) ){
		label_map <- aes_string(x=x, y=y, label=label, na.rm=TRUE)
		p <- p + geom_text(label_map, data=rm.na.phyloseq(DF, label),
					size=2, vjust=1.5, na.rm=TRUE)
	}

	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	# Return the ggplot object
	return(p)
}
################################################################################
# Define the ord.plot.DF.internal
################################################################################
#' @keywords internal
ord.plot.DF.internal <- function(physeq, ordination, type="sites", axes=c(1, 2)){

	coord <- scores(ordination, choices=axes, display=type)
	# coord row.names index order should match physeq. Enforce.
	if( type == "species" ){
		coord <- coord[taxa_names(physeq), ]
	} else if(type == "sites"){
		coord <- coord[sample_names(physeq), ]		
	}
	
	# If there is supplemental data, add it, else, return coord
	supp <- NULL
	# Define supplemental data. Use explicit accessor to avoid constructor options.
	if( !is.null(access(physeq, "sam_data")) & type == "sites"){
		supp  <- sample_data(physeq) # Supplemental data, samples
	} else if( !is.null(access(physeq, "tax_table")) & type == "species"){
		supp  <- tax_table(physeq) # Supplemental data, taxa
	}
	if( is.null(supp) ){
		DF <- coord
	} else {
		# Check that coord and supp have same indices. 
		if( !setequal(row.names(coord), row.names(supp)) ){
			stop("Ordination and supplementary data indices differ on the following:\n.",
				setdiff(row.names(coord), row.names(supp)))
		}
		# Combine for plotting data.frame
		DF <- data.frame(coord, supp)		
	}

	# Enforce DF class as data.frame.
	# Important in cases where no merging happens, scores may return a matrix, and then ggplot() fails.
	if( class(DF) != "data.frame"){ DF <- data.frame(DF) }
	
	return(DF)		
}
################################################################################
################################################################################
# Remove NA elements from data.frame prior to plotting
# Remove NA level from factor
################################################################################
#' @keywords internal
rm.na.phyloseq <- function(DF, key.var){
	# (1) Remove elements from DF if key.var has NA
	# DF[!is.na(DF[, key.var]), ]
	DF <- subset(DF, !is.na(eval(parse(text=key.var))))
	# (2) Remove NA from the factor level, if a factor.
	if( class(DF[, key.var]) == "factor" ){
		DF[, key.var] <- factor(as(DF[, key.var], "character"))
	}
	return(DF)
}
################################################################################
################################################################################
#' @keywords internal
#' @importFrom plyr is.discrete
rp.joint.fill <- function(DF, map.var, id.type.rp="samples"){
	# If all of the map.var values for samples/species are NA, replace with id.type.rp
	if( all(is.na(DF[DF$id.type==id.type.rp, map.var])) ){
		# If discrete, coerce to character, convert to factor, replace
		if( is.discrete(DF[, map.var]) ){
			temp.vec <- as(DF[, map.var], "character")
			temp.vec[is.na(temp.vec)] <- id.type.rp
			DF[, map.var] <- factor(temp.vec)
		}
	}
	return(DF)
}
################################################################################
################################################################################
#' Subset points from an ordination-derived ggplot
#'
#' Easily retrieve a plot-derived \code{data.frame} with a subset of points
#' according to a threshold and method. The meaning of the threshold depends
#' upon the method. See argument description below.
#' There are many useful examples of phyloseq ordination graphics in the
#' \href{http://joey711.github.com/phyloseq/subset_ord_plot-examples}{phyloseq online tutorials}.
#'
#' @usage subset_ord_plot(p, threshold=0.05, method="farthest")
#' 
#' @param p (Required).  A \code{\link{ggplot}} object created by 
#'  \code{\link{plot_ordination}}. It contains the complete data that you
#'  want to subset.
#'
#' @param threshold (Optional). A numeric scalar. Default is \code{0.05}.
#'  This value determines a coordinate threshold or population threshold,
#'  depending on the value of the \code{method} argument, ultimately 
#'  determining which points are included in returned \code{data.frame}.
#'
#' @param method (Optional). A character string. One of 
#'  \code{c("farthest", "radial", "square")}. Default is \code{"farthest"}.
#'  This determines how threshold will be interpreted.
#'
#' \describe{
#'
#'    \item{farthest}{
#'       Unlike the other two options, this option implies removing a 
#'       certain fraction or number of points from the plot, depending
#'       on the value of \code{threshold}. If \code{threshold} is greater
#'       than or equal to \code{1}, then all but \code{threshold} number 
#'       of points farthest from the origin are removed. Otherwise, if
#'       \code{threshold} is less than \code{1}, all but \code{threshold}
#'       fraction of points farthests from origin are retained.
#'    }
#' 
#'    \item{radial}{
#'	     Keep only those points that are beyond \code{threshold} 
#'       radial distance from the origin. Has the effect of removing a
#'       circle of points from the plot, centered at the origin.
#'    }
#' 
#'    \item{square}{
#'         Keep only those points with at least one coordinate
#'         greater than \code{threshold}. Has the effect of removing a 
#'         ``square'' of points from the plot, centered at the origin.
#'    }
#' 
#'  }
#'
#' @return A \code{\link{data.frame}} suitable for creating a 
#'  \code{\link{ggplot}} plot object, graphically summarizing
#'  the ordination result according to previously-specified parameters.
#' 
#' @seealso 
#'  \href{http://joey711.github.com/phyloseq/subset_ord_plot-examples}{phyloseq online tutorial} for this function.
#'
#'  \code{\link{plot_ordination}}
#'
#' @import ggplot2
#' @export
#' @examples 
#' ## See the online tutorials.
#' ## http://joey711.github.com/phyloseq/subset_ord_plot-examples
subset_ord_plot <- function(p, threshold=0.05, method="farthest"){
	threshold <- threshold[1] # ignore all but first threshold value.
	method    <- method[1] # ignore all but first string.
	method.names <- c("farthest", "radial", "square")
	# Subset to only some small fraction of points 
	# with furthest distance from origin
	df <- p$data[, c(1, 2)]
	d <- sqrt(df[, 1]^2 + df[, 2]^2)
	names(d) <- rownames(df)
	if( method.names[pmatch(method, method.names)] == "farthest"){
		if( threshold >= 1){
			show.names <- names(sort(d, TRUE)[1:threshold])
		} else if( threshold < 1 ){
			show.names <- names(sort(d, TRUE)[1:round(threshold*length(d))])
		} else {
			stop("threshold not a valid positive numeric scalar")
		}		
	} else if( method.names[pmatch(method, method.names)] == "radial"){
		show.names <- names(d[d > threshold])
	} else if( method.names[pmatch(method, method.names)] == "square"){	
		# show.names <- rownames(df)[as.logical((abs(df[, 1]) > threshold) + (abs(df[, 2]) > threshold))]
		show.names <- rownames(df)[((abs(df[, 1]) > threshold) | (abs(df[, 2]) > threshold))]
	} else {
		stop("method name not supported. Please select a valid method")
	}

	return(p$data[show.names, ])
}
################################################################################
################################################################################
#' General ordination eigenvalue plotter using ggplot2.
#'
#' Convenience wrapper for plotting ordination eigenvalues (if available) 
#' using a \code{ggplot2}-graphic.
#'
#' @param ordination (Required). An ordination object. Many different classes
#'  of ordination are defined by \code{R} packages. Ordination classes
#'  currently supported/created by the \code{\link{ordinate}} function are
#'  supported here.
#'  There is no default, as the expectation is that the 
#'  ordination will be performed and saved prior to calling this plot function.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @return A \code{\link{ggplot}} plot object, graphically summarizing
#'  the ordination result for the specified axes.
#' 
#' @seealso 
#'
#'  \code{\link{plot_ordination}}
#'
#'  \code{\link{ordinate}}
#'
#'  \code{\link{distance}}
#' 
#'  The examples on the phyloseq wiki page for \code{plot_ordination} show 
#'  many more examples:
#'
#' \url{https://github.com/joey711/phyloseq/wiki/plot_ordination}
#'
#' @import ggplot2
#' @export
#' @examples
#' # First load and trim a dataset
#' data(GlobalPatterns)
#' GP = prune_taxa(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
#' # Define a human-associated versus non-human categorical variable, and add new human variable to sample data:
#' sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )
#' # # filtering
#' # Remove taxa not seen more than 3 times in at least 20% of the samples
#' gp  = filter_taxa(GP, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
#' # Standardize abundances to the median sequencing depth
#' gpr = transform_sample_counts(gp, function(x, total=median(sample_sums(gp))) round(total * (x / sum(x))) )
#' # Let's use Coefficient of Variation for filtering, arbitrary cutoff of 3.0
#' gprf = filter_taxa(gpr, function(x) sd(x)/mean(x) > 3L, TRUE)
#' # For a somewhat readable number of taxa on display, let's subset to just Bacteroidetes for some plots
#' gprfb = subset_taxa(gprf, Phylum=="Bacteroidetes")
#' # Test plots (preforms ordination in-line, then makes scree plot)
#' plot_scree(ordinate(gprfb, "DPCoA", "bray"))
#' plot_scree(ordinate(gprfb, "PCoA", "bray"))
#' plot_scree(ordinate(gprfb, "NMDS", "bray")) # Empty return with message
#' plot_scree(ordinate(gprfb ~ SampleType, "CCA"))
#' plot_scree(ordinate(gprfb ~ SampleType, "RDA")) 
#' plot_scree(ordinate(gprfb, "DCA"))
#' plot_ordination(gprfb, ordinate(gprfb, "DCA"), type="scree")
plot_scree = function(ordination, title=NULL){
	# Use get_eigenvalue method dispatch. It always returns a numeric vector.
	x = extract_eigenvalue(ordination)
	# Were eigenvalues found? If not, return NULL
	if( is.null(x) ){
		cat("No eigenvalues found in ordination\n")
		return(NULL)
	} else {
		# If no names, add them arbitrarily "axis1, axis2, ..., axisN"
		if( is.null(names(x)) ) names(x) = 1:length(x)
		# For scree plot, want to show the fraction of total eigenvalues
		x = x/sum(x)
		# Set negative values to zero
		x[x <= 0.0] = 0.0		
		# Create the ggplot2 data.frame, and basic ggplot2 plot
		gdf = data.frame(axis=names(x), eigenvalue = x)
		p = ggplot(gdf, aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity")
		# Force the order to be same as original in x
		p = p + scale_x_discrete(limits = names(x))
		# Orient the x-labels for space.
		p = p + theme(axis.text.x = element_text(angle = 90))
		# Optionally add a title to the plot
		if( !is.null(title) ){
			p <- p + ggtitle(title)
		}		
		return(p)
	}
}
################################################################################
# Define S3 generic extract_eigenvalue function; formerly S4 generic get_eigenvalue()
# Function is used by `plot_scree` to get the eigenvalue vector from different
# types of ordination objects. 
# Used S3 generic in this case because many ordination objects, the input, are
# not formally-defined S4 classes, but vaguely-/un-defined S3. This throws
# warnings during package build if extract_eigenvalue were S4 generic method,
# because the ordination classes don't appear to have any definition in phyloseq
# or dependencies.
#' @keywords internal
extract_eigenvalue = function(ordination) UseMethod("extract_eigenvalue", ordination)
# Default is to return NULL (e.g. for NMDS, or non-supported ordinations/classes).
extract_eigenvalue.default = function(ordination) NULL
# for pcoa objects
extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
# for CCA objects
extract_eigenvalue.cca = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for RDA objects
extract_eigenvalue.rda = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for dpcoa objects
extract_eigenvalue.dpcoa = function(ordination) ordination$eig
# for decorana (dca) objects
extract_eigenvalue.decorana = function(ordination) ordination$evals
################################################################################
#' Melt phyloseq data object into large data.frame
#'
#' The psmelt function is a specialized melt function for melting phyloseq objects
#' (instances of the phyloseq class), usually for the purpose of graphics production
#' in ggplot2-based phyloseq-generated graphics. It relies heavily on the 
#' \code{\link[reshape]{melt}} and \code{\link{merge}} functions. Note that
#' ``melted'' phyloseq data is stored much less efficiently, and so RAM storage
#' issues could arise with a smaller dataset
#' (smaller number of samples/OTUs/variables) than one might otherwise expect.
#' For average-sized datasets, however, this should not be a problem.
#' Because the number of OTU entries has a large effect on the RAM requirement,
#' methods to reduce the number of separate OTU entries, for instance by
#' agglomerating based on phylogenetic distance using \code{\link{tipglom}},
#' can help alleviate RAM usage problems.
#' This function is made user-accessible for flexibility, but is also used 
#' extensively by plot functions in phyloseq.
#'
#' @usage psmelt(physeq)
#'
#' @param physeq (Required). An \code{\link{otu_table-class}} or 
#'  \code{\link{phyloseq-class}}. Function most useful for phyloseq-class.
#'
#' @return A \code{\link{data.frame}}-class table.
#'
#' @seealso
#'  \code{\link{plot_bar}}
#' 
#'  \code{\link[reshape]{melt}}
#'
#'  \code{\link{merge}}
#' 
#' @import reshape
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#' mdf = psmelt(gp.ch)
#' nrow(mdf)
#' ncol(mdf)
#' colnames(mdf)
#' head(rownames(mdf))
#' # Create a ggplot similar to
#' library("ggplot2")
#' p = ggplot(mdf, aes(x=SampleType, y=Abundance, fill=Genus))
#' p = p + geom_bar(color="black", stat="identity", position="stack")
#' print(p)
psmelt = function(physeq){
	
	# enforce orientation
	otutab = otu_table(physeq)
	if( !taxa_are_rows(otutab) ){
		otutab = t(otutab)	
	}
	mot <- as(otutab, "matrix")
	mdf <- melt(mot)
	colnames(mdf)[1] = "OTU"
	colnames(mdf)[2] = "Sample"
		
	# Merge the sample data.frame if present
	if( !is.null(sample_data(physeq, FALSE)) ){
		sdf = data.frame(sample_data(physeq))
		sdf$Sample = sample_names(physeq)
		# merge the sample-data and the melted otu table
		mdf = merge(mdf, sdf, by.x="Sample")
	}

	# Next merge taxonomy data
	if( !is.null(tax_table(physeq, FALSE)) ){
		tdf = data.frame(tax_table(physeq), OTU=taxa_names(physeq))
		mdf = merge(mdf, tdf, by.x="OTU")	
	}
	
	# Annotate the "value" column as the measured OTU "Abundance"
	colnames(mdf)[colnames(mdf)=="value"] = "Abundance"
	
	# Sort the entries by abundance
	mdf = mdf[order(mdf$Abundance, decreasing=TRUE), ]
		
	return(mdf)
}
################################################################################
################################################################################
#' A flexible, informative barplot phyloseq data
#'
#' There are many useful examples of phyloseq barplot graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
#' This function wraps \code{ggplot2} plotting, and returns a \code{ggplot2}
#'  graphic object
#' that can be saved or further modified with additional layers, options, etc.
#' The main purpose of this function is to quickly and easily create informative
#' summary graphics of the differences in taxa abundance between samples in
#' an experiment. 
#'
#' @usage plot_bar(physeq, x="Sample", y="Abundance", fill=NULL,
#'  title=NULL, facet_grid=NULL)
#'
#' @param physeq (Required). An \code{\link{otu_table-class}} or 
#'  \code{\link{phyloseq-class}}.
#'
#' @param x (Optional). Optional, but recommended, especially if your data
#'  is comprised of many samples. A character string.
#'  The variable in the melted-data that should be mapped to the x-axis.
#'  See \code{\link{psmelt}}, \code{\link{melt}},
#'  and \code{\link{ggplot}} for more details.
#' 
#' @param y (Optional). A character string.
#'  The variable in the melted-data that should be mapped to the y-axis.
#'  Typically this will be \code{"Abundance"}, in order to
#'  quantitatively display the abundance values for each OTU/group. 
#'  However, alternative variables could be used instead,
#'  producing a very different, though possibly still informative, plot.
#'  See \code{\link{psmelt}}, \code{\link{melt}},
#'  and \code{\link{ggplot}} for more details.
#'
#' @param fill (Optional). A character string. Indicates which sample variable
#'  should be used to map to the fill color of the bars. 
#'  The default is \code{NULL}, resulting in a gray fill for all bar segments.
#' 
#' @param facet_grid (Optional). A formula object.
#'  It should describe the faceting you want in exactly the same way as for 
#'  \code{\link[ggplot2]{facet_grid}}, 
#'  and is ulitmately provided to \code{\link{ggplot}}2 graphics.
#'  The default is: \code{NULL}, resulting in no faceting.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @return A \code{\link[ggplot2]{ggplot}}2 graphic object -- rendered in the graphical device
#'  as the default \code{\link[base]{print}}/\code{\link[methods]{show}} method.
#'
#' @seealso 
#'  \href{http://joey711.github.com/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
#'
#'  \code{\link{psmelt}}
#'
#'  \code{\link{ggplot}}
#' 
#'  \code{\link{qplot}}
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#' plot_bar(gp.ch)
#' plot_bar(gp.ch, fill="Genus")
#' plot_bar(gp.ch, x="SampleType", fill="Genus")
#' plot_bar(gp.ch, "SampleType", fill="Genus", facet_grid=~Family)
#' # See additional examples in the plot_bar online tutorial. Link above.
plot_bar = function(physeq, x="Sample", y="Abundance", fill=NULL,
	title=NULL, facet_grid=NULL){
		
	# Start by melting the data in the "standard" way using psmelt.
	mdf = psmelt(physeq)
	
	# Build the plot data structure
	p = ggplot(mdf, aes_string(x=x, y=y, fill=fill))

	# Add the bar geometric object. Creates a basic graphic. Basis for the rest.
	# Test weather additional
	p = p + geom_bar(stat="identity", position="stack", color="black")

	# By default, rotate the x-axis labels (they might be long)
	p = p + theme(axis.text.x=element_text(angle=-90, hjust=0))

	# Add faceting, if given
	if( !is.null(facet_grid) ){	
		p <- p + facet_grid(facet_grid)
	}
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	return(p)
}
################################################################################
################################################################################
# plot_tree section.  
# Includes core code borrowed with-permission from and attribution to the
# ggphylo package available (only) on GitHub
################################################################################
################################################################################
#' Extracts the parent node index for the given node. 
#'
#' Returns -1 if the node is root.
#' Return the index of the node directly above the given node.
#' Returns -1 if the given node is root.
#'
#' @param phylo, input phylo object
#' @param node, integer index of the node whose parent is desired
#' @return integer, the index of the parent node or -1 if the given node is root.
#' 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#'
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.parent.node <- function(phylo, node) {
  edge.index <- which(phylo$edge[,2]==node)
  node <- phylo$edge[edge.index,1]
  if (length(node)==0) {
    node <- -1
  }
  return(node)
}
################################################################################
#' Extracts the length of the branch above the given node.
#'
#' Returns 0 if the node is root.
#' 
#' @param phylo input phylo object
#' @param node integer, the node's index
#' @return numeric, the branch length of the edge leading to the given node.
#' May be NA.
#' 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.branch.length <- function(phylo, node) {
  edge.index <- which(phylo$edge[,2]==node)
  if (is.null(phylo$edge.length)) {
    return(NA)
  }
  bl <- phylo$edge.length[edge.index]
  if (length(bl)==0) {
    bl <- 0
  }
  return(bl)
}
################################################################################
#' Return a list of a node's children.
#'
#' Returns a list (not a vector!) of the node indices of the given
#' node's direct children. Returns (-1, -1) if the given node is a leaf.
#'
#' @param phylo, input phylo object
#' @param node, integer index of the node to test
#' @return list, a list containing the integer indices 
#' of the nodes directly beneath the given node.
#' 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.child.nodes <- function(phylo, node) {
  edge.indices <- which(phylo$edge[,1]==node)
  edge.indices <- sort(edge.indices)
  nodes <- phylo$edge[edge.indices,2]
  if (length(nodes)==0) {
    nodes <- list(c(-1,-1))
  } else {
    nodes <- list(nodes)
  }
  return(list(nodes))
}
################################################################################
#' @keywords internal
is.standard.layout <- function(x) {
  any(x %in% c('default', 'radial'))
}
################################################################################
#' Return length to root from node.
#'
#' Returns the length from the tree root to the given node. The input
#'  node can either be input as a node index or a node label.
#' 
#' @param phylo input phylo object
#' @param node integer or character. When integer, the node index; when character, the node label
#' @return numeric, the total branch length separating the tree root and the given node.
#' 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.length.to.root <- function(phylo, node) {
  tip.index <- node
  if (is.character(node)) {
    tip.index <- which(phylo$tip.label==node)
  }
  cur.node.b <- tip.index

  p.edges <- phylo$edge
  p.lengths <- phylo$edge.length

  if(is.null(p.lengths)) {
    p.lengths <- rep(1, length(p.edges[,1]))
  }

  length <- 0
  while(length(which(p.edges[,2]==cur.node.b)) > 0) {
    cur.edge.index <- which(p.edges[,2]==cur.node.b)
    cur.edge.length <- p.lengths[cur.edge.index]
    if (length(cur.edge.length) == 0 || is.na(cur.edge.length)) {
      cur.edge.length <- 0
    }
    length <- length + cur.edge.length
    cur.node.a <- p.edges[cur.edge.index,1]
    cur.node.b <- cur.node.a # Move up to the next edge
  }
  return(length)
}
################################################################################
#' Convert tree tags into data.frame columns
#' 
#' Given a \code{\link{phylo}} object and a data frame, transform all
#' the tags from the tree into columns of the data frame. Rows of the
#' data frame are linked to the tree via a required 'node' column, which
#' must contain integer indices of the associated node.
#'
#' This function is similar to the tree.as.data.frame method, but not
#' exactly the same. It is used internally by the tree.layout
#' function.
#' 
#' @param phylo, input phylo object
#' @param df, data.frame with a 'node' column corresponding to integer indices
#' of tree nodes.
#' @return df, a copy of the input data frame, with tags added as new columns
#'
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tags.into.df <- function(phylo, df) {
  all.tags <- c()

  for (i in 1:nrow(df)) {
    row <- df[i,]
    node <- row$node
    tags <- tree.get.tags(phylo, node)
    for (j in 1:length(tags)) {
      all.tags <- c(all.tags, names(tags)[j])
    }
  }
  all.tags <- unique(all.tags)

  df[, all.tags] <- 0

  for (i in 1:nrow(df)) {
    row <- df[i,]
    node <- row$node
    tags <- tree.get.tags(phylo, node)

    for (j in 1:length(tags)) {
      key <- names(tags)[j]
      val <- tags[j]
      df[i, key] <- val
    }
  }
  return(df)
}
################################################################################
#' Retrieves a list of all tags for the given node.
#'
#' @param phylo input phylo object
#' @param node the node index for the desired tags
#' @return list containing all tags associated with this node, if tags exist; empty list otherwise.
#'
#' @examples
#' # tree <- tree.read('((a,b[&&NHX:foo=bar]),c);')
#' # tree.get.tags(tree, tree.node.with.label(tree, 'b')) # foo => bar
#' # 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.get.tags <- function(phylo, node) {
  if (!tree.has.tags(phylo)) {
    return(list())
  }
  
  tags <- phylo$.tags[[node]]
  #print(paste(label(phylo, node), tags))
  if (is.null(tags) || is.na(tags)) {
    return(list())
  } else {
    return(tags)
  }
}
################################################################################
#' Determines whether the given phylo object contains tags or not.
#'
#' @param phylo input phylo object
#' @return boolean, indicating this phylo has tags (TRUE) or doesn't (FALSE).
#' 
#' @examples
#' # tree.has.tags(tree.read('((a,b[&&NHX:foo=bar]),c);')) # TRUE
#' # tree.has.tags(tree.read('((a,b),c);')) # FALSE
#' # 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @keywords internal
tree.has.tags <- function(phylo) {
  !is.null(phylo$.tags)
}
################################################################################
#' Returns a data frame defining segments to draw the phylogenetic tree.
#'
#' This internal function is borrowed directly from the \code{ggphylo} package
#' available on GitHub: \url{https://github.com/gjuggler/ggphylo}
#' 
#' @seealso
#' This code is borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' 
#' @author Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @importFrom plyr rbind.fill
#' @import ape
#'
#' @keywords internal
tree.layout <- function(
  phylo,
  layout = 'default',
  layout.ancestors = FALSE,
  ladderize=FALSE,
  align.seq.names = NA
) {

	if (ladderize != FALSE) {
		if (ladderize == 'left') {
			phylo <- ladderize(phylo, FALSE)
		} else {
			phylo <- ladderize(phylo, TRUE)
		}
	}

  # Number of nodes and leaves.
  n.nodes <- length(phylo$tip.label)+phylo$Nnode
  n.leaves <- length(phylo$tip.label)

  t.labels <- phylo$tip.label
  n.labels <- ((n.leaves+1):n.nodes)
  if (!is.null(phylo$node.label)) {
    n.labels <- phylo$node.label
  }

  # Create the skeleton data frame.
  # node     - Nodes with IDs 1 to N
  # x, y     - These will contain the x and y coords after the pending layout procedure.
  # label    - The first n.leaves nodes are the labeled tips
  # is.leaf  - Store is-leaf boolean for convenience
  # parent   - Contain the ID of the current node's parent
  # children - List of IDs of the current node's children
  # branch.length - Contains the branch lengths
  df <- data.frame(node     = c(1:n.nodes),
                   angle    = 0,
                   x        = 0,
                   y        = 0,
                   label    = c(t.labels, n.labels),
                   is.leaf  = c(rep(TRUE, n.leaves), rep(FALSE, n.nodes-n.leaves)),
                   parent   = 0,                                                     
                   children = 0,                                                  
                   branch.length = 0
		)

  # Collect the parents, children, and branch lengths for each node
  parent <- c()
  bl <- c()
  children <- list()
  event.count <- c()
  for (i in 1:nrow(df)) {
    node <- df[i,]$node
    parent <- c(parent, tree.parent.node(phylo, node))
    bl <- c(bl, tree.branch.length(phylo, node))
    children <- c(children, tree.child.nodes(phylo, node))
  }
  df$parent <- parent
  df$branch.length <- bl
  df$children <- children

  # Start the layout procedure by equally spacing the leaves in the y-dimension.
  # ape uses the edge ordering to indicate the plot position.
  # So we assign starting y-values according to the rank of each tip's location in
  # the edge vector.
  leaf.nodes <- which(df$is.leaf == TRUE)
  leaf.node.edge.indices <- match(leaf.nodes, phylo$edge[,2])
  df[df$is.leaf==TRUE,]$y <- rank(leaf.node.edge.indices)

  found.any.internal.node.sequences <- FALSE

  if (is.standard.layout(layout)) {
    # For each leaf: travel up towards the root, laying out each internal node along the way.
    for (i in 1:n.leaves) {
      cur.node <- i
      while (length(cur.node) > 0 && cur.node != -1) {
        df[cur.node, 'angle'] <- 0

        # We always use branch lengths: x-position is simply the length to the root.
        df[cur.node,]$x <- tree.length.to.root(phylo,cur.node)

        # The y-position for internal nodes is the mean of the y-position of all children.
        children <- unlist(df[cur.node,]$children)
        if (length(children) > 0 && children[1] != -1) {
          child.y.sum <- 0
          for (i in 1:length(children)) {
            child.index <- children[i]
            cur.y <- df[child.index,]$y
            child.y.sum <- child.y.sum + cur.y
          }
          df[cur.node, ]$y <- (child.y.sum) / length(children)
        }

        # Try to find the index of this node in the alignment names.
        if (!is.na(align.seq.names)) {
          lbl <- df[cur.node,]$label
          index.in.names <- which(align.seq.names == lbl | align.seq.names %in% c(paste('Node',lbl),
					          paste('Root node',lbl)))
          if (length(index.in.names)>0) {
            df[cur.node,]$y <- index.in.names
            if (!df[cur.node,]$is.leaf) {
              found.any.internal.node.sequences <- TRUE
            }
          }
        }

        cur.node <- unlist(df[cur.node,]$parent)
      }
    }
  }

	if (layout == 'unrooted') {
	# Not currently supported option.
	# 
    # # See http://code.google.com/p/phylowidget/source/browse/trunk/PhyloWidget/src/org/phylowidget/render/LayoutUnrooted.java
    # # For unrooted layout, we start from the root.
    # layout.f <- function(node, lo, hi) {
      # cur.enclosed <- tree.leaves.beneath(phylo, node)
      # cur.x <- df[node, 'x']
      # cur.y <- df[node, 'y']

      # children <- unlist(df[node, ]$children)
      # if (length(children) > 0 && children[1] != -1) {
        # cur.angle <- lo
        # for (i in 1:length(children)) {
          # child.node <- children[i]
          # child.enclosed <- tree.leaves.beneath(phylo, child.node)
          # child.ratio <- child.enclosed / cur.enclosed
          # child.bl <- tree.branch.length(phylo, child.node)

          # arc.size <- (hi - lo) * child.ratio
          # mid.angle <- cur.angle + arc.size / 2
          
          # child.x <- cur.x + cos(mid.angle) * child.bl
          # child.y <- cur.y + sin(mid.angle) * child.bl

          # df[child.node, 'x'] <<- child.x
          # df[child.node, 'y'] <<- child.y
          # df[child.node, 'angle'] <<- mid.angle / (2*pi) * 360

          # layout.f(child.node, cur.angle, cur.angle+arc.size)
          # cur.angle <- cur.angle + arc.size          
        # }        
      # }      
    # }
    
    # layout.f(tree.get.root(phylo), 0, 2 * pi)
  }

  df$dir <- 'none'

  # We have a data frame with each node positioned.
  # Now we go through and make two line segments for each node (for a 'square corner' type tree plot).
  line.df <- data.frame()
  for (i in 1:nrow(df)) {
    row <- df[i,]            # Data frame row for the current node.
    if (row$parent == -1) {
      next; # Root node!
    }
    p.row <- df[row$parent,] # Data frame row for the parent node.

    if (is.standard.layout(layout) && !(layout.ancestors && found.any.internal.node.sequences)) {
      horiz.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=row$y,
                               x=row$x,
                               xend=p.row$x,
                               label=row$label,   
                               dir='up',
                               branch.length=row$branch.length
                               )
      vert.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=p.row$y,
                               x=p.row$x,
                               xend=p.row$x,
                               label=row$label,
                               dir='across',
                               branch.length=row$branch.length
      )
      line.df <- rbind(line.df, horiz.line, vert.line)
    } else {
      up.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=p.row$y,
                               x=row$x,
                               xend=p.row$x,
                               label=row$label,
                               dir='up',
                               branch.length=row$branch.length
                               )
      line.df <- rbind(line.df, up.line)
    }
  }

  line.df <- tags.into.df(phylo, line.df)
  df <- tags.into.df(phylo, df)
  # Remove weird list-of-lists from the nodes data frame.
  df$children <- NULL

  label.df <- df
  line.df$type <- 'line'
  df$type <- 'node'
  label.df$type <- 'label'

  internal.label.df <- subset(label.df, is.leaf==FALSE)
  internal.label.df$type <- 'internal.label'

  label.df <- subset(label.df, is.leaf==TRUE)  

  all.df <- rbind.fill(line.df, df, label.df, internal.label.df)
  all.df
}
################################################################################
# Define an internal function for determining what the text-size should be
#' @keywords internal
manytextsize <- function(n, mins=0.5, maxs=4, B=6, D=100){
	# empirically selected size-value calculator.
	s <- B * exp(-n/D)
	# enforce a floor.
	s <- ifelse(s > mins, s, mins)
	# enforce a max
	s <- ifelse(s < maxs, s, maxs)
	return(s)
}
################################################################################
# Define an internal function for mapping phyloseq data variables to melted.tip
#' @keywords internal
treeMapVar2Tips <- function(melted.tip, physeq, variate){
	# If variate is tax_table-variable: Map tax_table-variable to melted.tip
	if( variate %in% rank_names(physeq, FALSE) ){
		# Add relevant tax_table column.
		x <- as(tax_table(physeq), "matrix")[, variate, drop=TRUE]
		names(x) <- taxa_names(physeq)
		return( x[as(melted.tip$taxa_names, "character")] )
	}
	# If variate is sampleMap-variable: Map sample-variable to melted.tip
	if( variate %in% sample_variables(physeq, FALSE) ){
		x <- as.vector(data.frame(sample_data(physeq))[, variate])
		names(x) <- sample_names(physeq)				
		return( x[as(melted.tip$variable, "character")] )
	}	
}
################################################################################
# The "tree only" setting. Simple. No annotations.
#' @keywords internal
#' @import ggplot2
plot_tree_only <- function(tdf){
	# build tree lines
	p <- ggplot(subset(tdf, type == "line")) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
	# Return ggplot object
	return(p)
}
################################################################################
# The "sampledodge" plot_tree subset function.
# Assumes the tree data.frame, tdf, has already been built and is third argument.
#' @keywords internal
#' @import ggplot2
#' @import reshape 
#' @import scales
#' @importFrom plyr aaply
#' @importFrom plyr ddply
plot_tree_sampledodge <- function(physeq, p, tdf, color, shape, size, min.abundance, 
				label.tips, text.size, sizebase, base.spacing){
								
	# Get the subset of tdf for just the tips (leaves)
	speciesDF <- subset(tdf, type=="label")
	
	# Add abundance data for each species
	# # First, re-order speciesDF ensure match with otu_table
	rownames(speciesDF) <- as(speciesDF$label, "character")
	speciesDF <- speciesDF[taxa_names(physeq), ]
	# # subset speciesDF to just what you need for tip plotting
	speciesDF <- data.frame(speciesDF[, c("x", "y")], taxa_names=rownames(speciesDF))
	
	# # Make the 0-values NA so they're not plotted. 
	OTU 		<- as(otu_table(physeq), "matrix") # Coerce to matrix.
	if(!taxa_are_rows(physeq)){OTU <- t(OTU)} # Enforce orientation.
	OTU[OTU==0] <- NA
	
	# # Now add abundance table
	speciesDF 	<- cbind(speciesDF, OTU)
	
	# # Now melt to just what you need for adding to plot
	melted.tip <- melt.data.frame(speciesDF, id.vars=c("x", "y", "taxa_names"))
	
	# Determine the horizontal adjustment index for each point
	h.adj <- aaply(OTU, 1, function(j){ 1:length(j) - cumsum(is.na(j)) - 1 })
	# Add to melted data.frame (melted.tip) - uses melt.array
	melted.tip$h.adj.index <- melt(h.adj)$value
	
	# the amount to adjust the horizontal coordinate space
	x.spacer.base    <- base.spacing * max(melted.tip$x)
	melted.tip$x.spacer.base <- x.spacer.base
	melted.tip$x.adj <- (melted.tip$h.adj.index * melted.tip$x.spacer.base)
	
	# Remove the NA values (the samples that had no individuals of a particular species)
	melted.tip <- subset(melted.tip, !is.na(value))
	if( nrow(melted.tip)==0L ){
		stop("The number of rows of tip data.frame has dropped to 0 after rm NA values")		
	}

	# Build the tip-label portion of the melted.tip data.frame, if needed.
	if( !is.null(label.tips) ){
		if( label.tips == "taxa_names" ){
			melted.tip$tipLabels <- melted.tip[, "taxa_names"]
		} else {
			melted.tip$tipLabels <- treeMapVar2Tips(melted.tip, physeq, label.tips)
		}
	} 	

	# color-map handling. Names "variable", "value" have specieal meaning.	
	if( !is.null(color) ){
		if( color %in% c("sample_names", "samples") ){
			color <- "variable"
		} else {
			melted.tip$color <- treeMapVar2Tips(melted.tip, physeq, color)
			names(melted.tip)[names(melted.tip)=="color"] <- color # rename to name of color variable
		}
	}

	# shape-map handling. Names "variable", "value" have specieal meaning.	
	if( !is.null(shape) ){
		if( shape %in% c("sample_names", "samples") ){
			shape <- "variable"
		} else if( !is.null(shape) ){
			melted.tip$shape <- treeMapVar2Tips(melted.tip, physeq, shape)
			names(melted.tip)[names(melted.tip)=="shape"] <- shape # rename to name of shape variable
		}
	}
	
	# size-map handling. Names "abundance", "variable", "value" have special meaning.
	ab_labels = c("abundance", "Abundance", "abund")
	if( !is.null(size) ){	
		if( size %in% ab_labels ){
			size = "value"
		} else {
			melted.tip$size <- treeMapVar2Tips(melted.tip, physeq, size)
			names(melted.tip)[names(melted.tip)=="size"] <- size # rename to name of size variable
		}
	}
		
	# The general tip-point map. Objects can be NULL, and that aesthetic gets ignored.
	tip.map <- aes_string(x="x + x.adj + x.spacer.base", y="y", color=color, fill=color, shape=shape, size=size)
	
	# Add the new point layer.
	p <- p + geom_point(tip.map, data=melted.tip, na.rm=TRUE)

	# Optionally-add abundance value label to each point.
	# This size needs to match point size.
	if( any(melted.tip$value >= min.abundance[1]) ){
		if( is.null(size) ){
			point.label.map <- aes_string(x="x + x.adj + x.spacer.base", y="y", label="value")
			p <- p + geom_text( point.label.map, data=subset(melted.tip, value>=min.abundance[1]), size=1, na.rm=TRUE)
		} else {
			point.label.map <- aes_string(x="x + x.adj + x.spacer.base", y="y",
				label="value", size=paste("0.025*", size, sep=""))
			p <- p + geom_text( point.label.map, angle=45, hjust=0,
						data=subset(melted.tip, value>=min.abundance[1]), na.rm=TRUE)
		}
	}

	# If indicated, add the species labels to the right of points.
	if( !is.null(label.tips) ){
		# melted.tip.far has only one row per tip,
		# the farthest horiz. adjusted position (one for each taxa)
		melted.tip.far <- ddply(melted.tip, "taxa_names", function(df){
			df[df$h.adj.index == max(df$h.adj.index), , drop=FALSE]
		})
		# Create the tip-label aesthetic map.
		label.map <- aes(x=x + x.adj + 2*x.spacer.base, y=y, label=tipLabels)
		# Add labels layer to plotting object.
		p <- p + geom_text(label.map, data=melted.tip.far, size=I(text.size), hjust=0, na.rm=TRUE)
	}
	
	# Adjust point size transform
	if( !is.null(size) ){
		p <- p + scale_size_continuous(trans=log_trans(sizebase))
	}
	
	# Update legend-name of color or shape or size
	if( identical(color, "variable") ){
		p <- update_labels(p, list(colour = "Samples"))
	}
	if( identical(shape, "variable") ){
		p <- update_labels(p, list(shape  = "Samples"))
	}
	if( identical(size, "value") ){
		p <- update_labels(p, list(size = "Abundance"))
	}
			
	return(p)		
}
################################################################################
# Return TRUE if the nodes of the tree in the phyloseq object provided are unlabeled.
#' @keywords internal
nodesnotlabeled = function(physeq){
	if( is.null(phy_tree(physeq, FALSE)) ){
		warning("There is no phylogenetic tree in the object you have provided. Try phy_tree(physeq) to see.")
		return(TRUE)
	} else {
		return( is.null(phy_tree(physeq)$node.label) | length(phy_tree(physeq)$node.label)==0L )
	}
}
# A quick test function to decide how nodes should be labeled by default, if at all.
#  
#' @keywords internal
howtolabnodes = function(physeq){
	if(!nodesnotlabeled(physeq)){
		return(nodeplotdefault(manytextsize(ntaxa(physeq))))
	} else {
		return(nodeplotblank)
	}
}
################################################################################
#' Function to avoid plotting node labels
#'
#' Unlike, \code{\link{nodeplotdefault}} and \code{\link{nodeplotboot}},
#' this function does not return a function, but instead is provided
#' directly to the \code{nodelabf} argument of \code{\link{plot_tree}} to 
#' ensure that node labels are not added to the graphic.
#' Please note that you do not need to create or obtain the arguments to 
#' this function. Instead, you can provide this function directly to 
#' \code{\link{plot_tree}} and it will know what to do with it. Namely,
#' use it to avoid plotting any node labels.
#'
#' @usage nodeplotblank(p, nodelabdf)
#'
#' @param p (Required). The \code{\link{plot_tree}} graphic.
#'
#' @param nodelabdf (Required). The \code{data.frame} produced internally in 
#' \code{link{plot_tree}} to use as data for creating ggplot2-based tree graphics.  
#'
#' @return The same input object, \code{p}, provided as input. Unmodified.
#'
#' @seealso 
#' \code{\link{nodeplotdefault}}
#'
#' \code{\link{nodeplotboot}}
#'
#' \code{\link{plot_tree}}
#'
#' @import ggplot2
#' @export
#' @examples
#' data("esophagus")
#' plot_tree(esophagus)
#' plot_tree(esophagus, nodelabf=nodeplotblank)
nodeplotblank = function(p, nodelabdf){
	return(p)
}
################################################################################
#' Generates a function for labeling bootstrap values on a phylogenetic tree.
#'
#' Is not a labeling function itself, but returns one.
#' The returned function is specialized for labeling bootstrap values.
#' Note that the function that 
#' is returned has two completely different arguments from the four listed here:
#' the plot object already built by earlier steps in
#' \code{\link{plot_tree}}, and the \code{\link{data.frame}}
#' that contains the relevant plotting data for the nodes
#' (especially \code{x, y, label}),
#' respectively.  
#' See \code{\link{nodeplotdefault}} for a simpler example.
#' The main purpose of this and \code{\link{nodeplotdefault}} is to
#' provide a useful default function generator for arbitrary and
#' bootstrap node labels, respectively, and also to act as 
#' examples of functions that can successfully interact with 
#' \code{\link{plot_tree}} to add node labels to the graphic.
#'
#' @usage nodeplotboot(highthresh=95L, lowcthresh=50L, size=2L, hjust=-0.2)
#'
#' @param highthresh (Optional). A single integer between 0 and 100.
#'  Any bootstrap values above this threshold will be annotated as
#'  a black filled circle on the node, rather than the bootstrap
#'  percentage value itself.
#'
#' @param lowcthresh (Optional). A single integer between 0 and 100,
#'  less than \code{highthresh}. Any bootstrap values below this value
#'  will not be added to the graphic. Set to 0 or below to add all
#'  available values.
#'
#' @param size (Optional). Numeric. Should be positive. The 
#'  size parameter used to control the text size of taxa labels.
#'  Default is \code{2}. These are ggplot2 sizes.
#'
#' @param hjust (Optional). The horizontal justification of the
#'  node labels. Default is \code{-0.2}.  
#'
#' @return A function that can add a bootstrap-values layer to the tree graphic.
#'  The values are represented in two ways; either as black filled circles
#'  indicating very high-confidence nodes, or the bootstrap value itself
#'  printed in small text next to the node on the tree.
#'
#' @seealso 
#' \code{\link{nodeplotdefault}}
#'
#' \code{\link{nodeplotblank}}
#'
#' \code{\link{plot_tree}}
#'
#' @import ggplot2
#' @export
#' @examples
#' nodeplotboot()
#' nodeplotboot(3, -0.4)
nodeplotboot = function(highthresh=95L, lowcthresh=50L, size=2L, hjust=-0.2){
	function(p, nodelabdf){
		# For bootstrap, check that the node labels can be coerced to numeric
		try(boot <- as(as(nodelabdf$label, "character"), "numeric"), TRUE)
		# Want NAs/NaN to propagate, but still need to test remainder
		goodboot = boot[complete.cases(boot)]
		if( !is(goodboot, "numeric") & length(goodboot) > 0 ){
			stop("The node labels, phy_tree(physeq)$node.label, are not coercable to a numeric vector with any elements.")
		}
		# So they look even more like bootstraps and display well, 
		# force them to be between 0 and 100, rounded to 2 digits.
		if( all( goodboot >= 0.0 & goodboot <= 1.0 ) ){
			boot = round(boot, 2)*100L
		}
		nodelabdf$boot = boot
		boottop = subset(nodelabdf, boot >= highthresh)
		bootmid = subset(nodelabdf, boot > lowcthresh & boot < highthresh)
		# Label the high-confidence nodes with a point.
		if( nrow(boottop)>0L ){
			p = p + geom_point(data = boottop, aes(x=x, y=y), na.rm=TRUE)
		}
		# Label the remaining bootstrap values as text at the nodes.
		if( nrow(bootmid)>0L ){
			bootmid$label = bootmid$boot
			p = nodeplotdefault(size, hjust)(p, bootmid)
		}
		return(p)
	}
}
################################################################################
#' Generates a default node-label function 
#'
#' Is not a labeling function itself, but returns one.
#' The returned function is capable of adding
#' whatever label is on a node. Note that the function that 
#' is returned has two completely different arguments to those listed here:
#' the plot object already built by earlier steps in
#' \code{\link{plot_tree}}, and the \code{\link{data.frame}}
#' that contains the relevant plotting data for the nodes
#' (especially \code{x, y, label}),
#' respectively. 
#' See \code{\link{nodeplotboot}} for a more sophisticated example.
#' The main purpose of this and \code{\link{nodeplotboot}} is to
#' provide a useful default function generator for arbitrary and
#' bootstrap node labels, respectively, and also to act as 
#' examples of functions that will successfully interact with 
#' \code{\link{plot_tree}} to add node labels to the graphic.
#'
#' @usage nodeplotdefault(size=2L, hjust=-0.2)
#'
#' @param size (Optional). Numeric. Should be positive. The 
#'  size parameter used to control the text size of taxa labels.
#'  Default is \code{2}. These are ggplot2 sizes.
#'
#' @param hjust (Optional). The horizontal justification of the
#'  node labels. Default is \code{-0.2}.  
#'
#' @return A function that can add a node-label layer to a graphic.
#'
#' @seealso 
#' \code{\link{nodeplotboot}}
#'
#' \code{\link{nodeplotblank}}
#'
#' \code{\link{plot_tree}}
#'
#' @import ggplot2
#' @export
#' @examples
#' nodeplotdefault()
#' nodeplotdefault(3, -0.4)
nodeplotdefault = function(size=2L, hjust=-0.2){
	function(p, nodelabdf){
		p = p + geom_text(data=nodelabdf, aes(x=x, y=y, label=label), size=size, hjust=hjust, na.rm=TRUE)
		return(p)
	}
}
################################################################################
#' Plot a phylogenetic tree with optional annotations
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
#' This function is intended to facilitate easy graphical investigation of 
#' the phylogenetic tree, as well as sample data. Note that for phylogenetic
#' sequencing of samples with large richness, some of the options in this 
#' function will be prohibitively slow to render, or too dense to be
#' interpretable. A rough ``rule of thumb'' is to use subsets of data 
#' with not many more than 200 OTUs per plot, sometimes less depending on the
#' complexity of the additional annotations being mapped to the tree. It is 
#' usually possible to create an unreadable, uninterpretable tree with modern
#' datasets. However, the goal should be toward parameter settings and data
#' subsets that convey (honestly, accurately) some biologically relevant
#' feature of the data. One of the goals of the \code{\link{phyloseq-package}}
#' is to make the determination of these features/settings as easy as possible.
#'
#' This function received a development contribution from the work of 
#' Gregory Jordan for the \code{ggphylo} package, which provides tools for 
#' rendering a phylogenetic tree in \code{\link{ggplot2}} graphics. That package
#' is not currently available from CRAN or Bioconductor, but is available 
#' in development (roughly ``beta'') form from GitHub. 
#' Furthermore, although \code{ggphylo} awesomely supports radial and unrooted trees, 
#' \code{plot_tree} currently only supports ``standard'' square-horizontal trees.
#' Additional support for these types of features (like radial trees)
#' is planned. Send us development feedback if this is a feature you really
#' want to have soon.
#'
#' @usage plot_tree(physeq, method="sampledodge", nodelabf=NULL, color=NULL, shape=NULL, size=NULL,
#'  min.abundance=Inf, label.tips=NULL, text.size=NULL, sizebase=5, base.spacing=0.02,
#' 	ladderize=FALSE, plot.margin=0.2, title=NULL)
#'
#' @param physeq (Required). The data about which you want to 
#'  plot and annotate a phylogenetic tree, in the form of a
#'  single instance of the \code{\link{phyloseq-class}}, containing at 
#'  minimum a phylogenetic tree component (try \code{\link{phy_tree}}).
#'  One of the major advantages of this function over basic tree-plotting utilities
#'  in the \code{\link{ape}}-package is the ability to easily annotate the tree
#'  with sample variables and taxonomic information. For these uses, 
#'  the \code{physeq} argument should also have a \code{\link{sample_data}}
#'  and/or \code{\link{tax_table}} component(s).
#' 
#' @param method (Optional). Character string. Default \code{"sampledodge"}. 
#'  The name of the annotation method to use. 
#'  This will be expanded in future versions.
#'  Currently only \code{"sampledodge"} and \code{"treeonly"} are supported.
#'  The \code{"sampledodge"} option results in points
#'  drawn next to leaves if individuals from that taxa were observed,
#'  and a separate point is drawn for each sample.
#' 
#' @param nodelabf (Optional). A function. Default \code{NULL}.
#'  If \code{NULL}, the default, a function will be selected for you based upon
#'  whether or not there are node labels in \code{phy_tree(physeq)}.
#'  For convenience, the phyloseq package includes two generator functions
#'  for adding arbitrary node labels (can be any character string),
#'  \code{\link{nodeplotdefault}};
#'  as well as for adding bootstrap values in a certain range,
#'  \code{\link{nodeplotboot}}.
#'  To not have any node labels in the graphic, set this argument to
#'  \code{\link{nodeplotblank}}.
#'
#' @param color (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point color.
#' 
#' @param shape (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point shape.
#'
#' @param size (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point size.
#'  A special argument \code{"abundance"} is reserved here and scales
#'  point size using abundance in each sample on a log scale.
#'
#' @param min.abundance (Optional). Numeric. 
#'  The minimum number of individuals required to label a point
#'  with the precise number.
#'  Default is \code{Inf},
#'  meaning that no points will have their abundance labeled.
#'  If a vector, only the first element is used. 
#'
#' @param label.tips (Optional). Character string. Default is \code{NULL},
#'  indicating that no tip labels will be printed.
#'  If \code{"taxa_names"}, then the name of the taxa will be added 
#'  to the tree; either next to the leaves, or next to
#'  the set of points that label the leaves. Alternatively,
#'  if this is one of the rank names (from \code{rank_names(physeq)}),
#'  then the identity (if any) for that particular taxonomic rank
#'  is printed instead.
#'
#' @param text.size (Optional). Numeric. Should be positive. The 
#'  size parameter used to control the text size of taxa labels.
#'  Default is \code{NULL}. If left \code{NULL}, this function
#'  will automatically calculate a (hopefully) optimal text size
#'  given the vertical constraints posed by the tree itself. 
#'  This argument is included in case the 
#'  automatically-calculated size is wrong, and you want to change it.
#'  Note that this parameter is only meaningful if \code{label.tips}
#'  is not \code{NULL}.
#'
#' @param sizebase (Optional). Numeric. Should be positive.
#'  The base of the logarithm used
#'  to scale point sizes to graphically represent abundance of
#'  species in a given sample. Default is 5.
#' 
#' @param base.spacing (Optional). Numeric. Default is \code{0.02}.
#'  Should be positive.
#'  This defines the base-spacing between points at each tip/leaf in the
#'  the tree. The larger this value, the larger the spacing between points.
#'  This is useful if you have problems with overlapping large points
#'  and/or text indicating abundance, for example. Similarly, if you 
#'  don't have this problem and want tighter point-spacing, you can 
#'  shrink this value.
#'
#' @param ladderize (Optional). Boolean or character string (either
#'  \code{FALSE}, \code{TRUE}, or \code{"left"}). Default is \code{FALSE}.
#'  This parameter specifies whether or not to \code{\link[ape]{ladderize}} the tree 
#'  (i.e., reorder nodes according to the depth of their enclosed
#'  subtrees) prior to plotting. When set to \code{TRUE}, the default
#'  ladderization (``right'' ladderization) is used; when set to
#'  \code{FALSE}, no ladderization is performed; when set to \code{"left"},
#'  the reverse direction (``left'' ladderization) is applied.
#'
#' @param plot.margin (Optional). Numeric. Default is \code{0.2}.
#'  Should be positive.
#'  This defines how much right-hand padding to add to the tree plot,
#'  which can be required to not truncate tip labels. The margin value
#'  is specified as a fraction of the overall tree width which is added
#'  to the right side of the plot area. So a value of \code{0.2} adds
#'  twenty percent extra space to the right-hand side of the plot.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @return A \code{\link{ggplot}}2 plot.
#' 
#' @seealso
#'  \code{\link{plot.phylo}}
#'
#' This function is a special use-case that relies 
#' on code borrowed directly, with permission, from the
#' not-yet-officially-released package, \code{ggphylo}, currently only
#' available from GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
#'
#' @author Paul McMurdie, relying on supporting code from
#'  Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @import reshape
#' @import scales
#' @import ggplot2
#' @export
#' @examples
#' # # Using plot_tree() with the esophagus dataset.
#' # # Please note that many more interesting examples are shown
#' # # in the online tutorials"
#' # # http://joey711.github.com/phyloseq/plot_tree-examples
#' data(esophagus)
#' plot_tree(esophagus)
#' plot_tree(esophagus, color="samples")
#' plot_tree(esophagus, size="abundance")
#' plot_tree(esophagus, size="abundance", color="samples")
#' plot_tree(esophagus, size="abundance", color="samples", base.spacing=0.03)
plot_tree <- function(physeq, method="sampledodge", nodelabf=NULL,
	color=NULL, shape=NULL, size=NULL,
	min.abundance=Inf, label.tips=NULL, text.size=NULL,
	sizebase=5, base.spacing = 0.02,
	ladderize=FALSE, plot.margin=0.2, title=NULL){

	# Test that physeq has tree, top-level test.
	if( is.null(phy_tree(physeq, FALSE)) ){
		stop("There is no phylogenetic tree in the object you have provided. Try phy_tree(physeq) to see.")
	}

	# Create the tree data.frame
	tdf <- tree.layout(phy_tree(physeq), ladderize=ladderize)

	# "Naked" unannotated tree built in ggplot2 no matter what. Lines only.
	p <- plot_tree_only(tdf)
	
	# If no text.size given, calculate it from number of tips ("species", aka taxa)
	# This is very fast. No need to worry about whether text is printed or not. DRY.
	if( is.null(text.size) ){
		text.size <- manytextsize(ntaxa(physeq))
	}

	# Tip annotation section.
	#
	# Annotate dodged sample points, and other fancy tip labels
	if( method == "sampledodge" ){
		p <- plot_tree_sampledodge(physeq, p, tdf, color, shape, size, min.abundance, 
				label.tips, text.size, sizebase, base.spacing)
	}

	# Node label section.
	# 
	# If no nodelabf ("node label function") given, ask internal function to pick one.
	# Is NULL by default, meaning will dispatch to howtolabnodes to select function.
	# For no node labels, the "dummy" function nodeplotnot will return tree plot 
	# object, p, as-is, unmodified.
	if( is.null(nodelabf) ){
		nodelabf = howtolabnodes(physeq)
	}
	# Subset data.frame to just the internal nodes (not leaves).
	nodelabdf = subset(tdf, is.leaf == FALSE & type == "node")
	# Use the provided/inferred node label function to add the node labels layer(s)
	p = nodelabf(p, nodelabdf)
	
	# Plot margins. 
	#
	# Adjust the tree graphic plot margins.
	# Helps to manually ensure that graphic elements aren't clipped,
	# especially when there are long tip labels.
	if( method == "sampledodge" ){
		min.x <- min(tdf$x, p$layers[[2]]$data$x, na.rm=TRUE)
		max.x <- max(tdf$x, p$layers[[2]]$data$x, na.rm=TRUE)
	} else {
		min.x <- min(tdf$x, na.rm=TRUE)
		max.x <- max(tdf$x, na.rm=TRUE)
	}
	if (plot.margin > 0) {
		max.x <- max.x * (1.0 + plot.margin)
	} 
	p <- p + scale_x_continuous(limits=c(min.x, max.x))	
	
	# Themeing section.
	#
	# Theme-ing: Blank theming 
	# Should open this up as function-argument also.
	p <- p + theme(axis.ticks = element_blank(),
			axis.title.x=element_blank(), axis.text.x=element_blank(),
			axis.title.y=element_blank(), axis.text.y=element_blank(),
			panel.background = element_blank(),
			panel.grid.minor = element_blank(),			
			panel.grid.major = element_blank()
			)
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	return(p)
}
################################################################################
################################################################################
################################################################################
# Adapted from NeatMap-package.
# Vectorized for speed and simplicity, also only calculates theta and not r.
#' @keywords internal
RadialTheta <- function(pos){
    pos = as(pos, "matrix")
    xc  = mean(pos[, 1])
    yc  = mean(pos[, 2])
    theta = atan2((pos[, 2] - yc), (pos[, 1] - xc))
    names(theta) <- rownames(pos)
    return(theta)
}
################################################################################
#' Create an ecologically-organized heatmap using ggplot2 graphics
#'
#' There are many useful examples of phyloseq heatmap graphics in the
#' \href{http://joey711.github.com/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
#' In a 2010 article in BMC Genomics, Rajaram and Oono show describe an 
#' approach to creating a heatmap using ordination methods to organize the 
#' rows and columns instead of (hierarchical) cluster analysis. In many cases
#' the ordination-based ordering does a much better job than h-clustering. 
#' An immediately useful example of their approach is provided in the NeatMap
#' package for R. The NeatMap package can be used directly on the abundance 
#' table (\code{\link{otu_table-class}}) of phylogenetic-sequencing data, but 
#' the NMDS or PCA ordination options that it supports are not based on ecological
#' distances. To fill this void, phyloseq provides the \code{plot_heatmap()}
#' function as an ecology-oriented variant of the NeatMap approach to organizing
#' a heatmap and build it using ggplot2 graphics tools.
#' The \code{distance} and \code{method} arguments are the same as for the
#' \code{\link{plot_ordination}} function, and support large number of
#' distances and ordination methods, respectively, with a strong leaning toward
#' ecology.
#' This function also provides the options to re-label the OTU and sample 
#' axis-ticks with a taxonomic name and/or sample variable, respectively, 
#' in the hope that this might hasten your interpretation of the patterns
#' (See the \code{sample.label} and \code{taxa.label} documentation, below). 
#' Note that this function makes no attempt to overlay hierarchical 
#' clustering trees on the axes, as hierarchical clustering is not used to 
#' organize the plot. Also note that each re-ordered axis repeats at the edge,
#' and so apparent clusters at the far right/left or top/bottom of the 
#' heat-map may actually be the same. For now, the placement of this edge
#' can be considered arbitrary, so beware of this artifact of this graphical
#' representation. If you benefit from this phyloseq-specific implementation
#' of the NeatMap approach, please cite both our packages/articles.
#'
#' This approach borrows heavily from the \code{heatmap1} function in the
#' \code{NeatMap} package. Highly recommended, and we are grateful for their
#' package and ideas, which we have adapted for our specific purposes here,
#' but did not use an explicit dependency. At the time of the first version
#' of this implementation, the NeatMap package depends on the rgl-package,
#' which is not needed in phyloseq, at present. Although likely a transient
#' issue, the rgl-package has some known installation issues that have further
#' influenced to avoid making NeatMap a formal dependency (Although we love
#' both NeatMap and rgl!).
#' 
#' @usage plot_heatmap(physeq, method="NMDS", distance="bray", 
#'  sample.label=NULL, taxa.label=NULL,
#'  low="#000033", high="#66CCFF", na.value="black", trans=log_trans(4),
#'  max.label=250, title=NULL, species.label=taxa.label, ...)
#'
#' @param physeq (Required). The data, in the form of an instance of the
#'  \code{\link{phyloseq-class}}. This should be what you get as a result
#'  from one of the
#'  \code{\link{import}} functions, or any of the processing downstream.
#'  No data components beyond the \code{\link{otu_table}} are strictly 
#'  necessary, though they may be useful if you want to re-label the 
#'  axis ticks according to some observable or taxonomic rank, for instance,
#'  or if you want to use a \code{\link{UniFrac}}-based distance
#'  (in which case your \code{physeq} data would need to have a tree included).
#'
#' @param method (Optional).
#'  The ordination method to use for organizing the 
#'  heatmap. A great deal of the usefulness of a heatmap graphic depends upon 
#'  the way in which the rows and columns are ordered. 
#'
#' @param distance (Optional). A character string. 
#'  The ecological distance method to use in the ordination.
#'  See \code{\link{distance}}.
#'
#' @param sample.label (Optional). A character string.
#'  The sample variable by which you want to re-label the sample (horizontal) axis.
#'
#' @param taxa.label (Optional). A character string.
#'  The name of the taxonomic rank by which you want to
#'  re-label the taxa/species/OTU (vertical) axis.
#'  You can see available options in your data using
#'  \code{\link{rank_names}(physeq)}.
#'
#' @param low (Optional). A character string. An R color.
#'  See \code{?\link{colors}} for options support in R (there are lots).
#'  The color that represents the lowest non-zero value
#'  in the heatmap. Default is a dark blue color, \code{"#000033"}.
#' 
#' @param high (Optional). A character string. An R color.
#'  See \code{\link{colors}} for options support in R (there are lots).
#'  The color that will represent the highest 
#'  value in the heatmap. The default is \code{"#66CCFF"}.
#'  Zero-values are treated as \code{NA}, and set to \code{"black"}, to represent
#'  a background color.
#'
#' @param na.value (Optional). A character string. An R color.
#'  See \code{\link{colors}} for options support in R (there are lots).
#'  The color to represent what is essentially the background of the plot,
#'  the non-observations that occur as \code{NA} or
#'  \code{0} values in the abundance table. The default is \code{"black"}, which 
#'  works well on computer-screen graphics devices, but may be a poor choice for
#'  printers, in which case you might want this value to be \code{"white"}, and
#'  reverse the values of \code{high} and \code{low}, above.
#'
#' @param trans (Optional). \code{"trans"}-class transformer-definition object.
#'  A numerical transformer to use in 
#'  the continuous color scale. See \code{\link[scales]{trans_new}} for details.
#'  The default is \code{\link{log_trans}(4)}.
#'
#' @param max.label (Optional). Integer. Default is 250.
#'  The maximum number of labeles to fit on a given axis (either x or y). 
#'  If number of taxa or samples exceeds this value, 
#'  the corresponding axis will be stripped of any labels. 
#'
#'  This supercedes any arguments provided to
#'  \code{sample.label} or \code{taxa.label}.
#'  Make sure to increase this value if, for example,
#'  you want a special label
#'  for an axis that has 300 indices.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'  
#' @param species.label (Deprecated). Equivalent to and over-ridden by
#'  \code{taxa.label}, but for the same purpose. Old nomenclature that
#'  will be removed in the next release of phyloseq. Included here for
#'  backward compatibility.
#'
#' @param ... (Optional). Additional parameters passed to \code{\link{ordinate}}.
#' 
#' @return
#'  A heatmap plot, in the form of a \code{\link{ggplot}2} plot object,
#'  which can be further saved and modified.
#'
#' @references
#'  Because this function relies so heavily in principle, and in code, on some of the
#'  functionality in NeatMap, please site their article if you use this function
#'  in your work.
#' 
#'  Rajaram, S., & Oono, Y. (2010).
#'  NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' Please see further examples in the 
#' \href{http://joey711.github.com/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
#'
#' @import vegan
#' @import scales 
#' 
#' @export
#' @examples
#' data("GlobalPatterns")
#' gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
#' # FYI, the base-R function uses a non-ecological ordering scheme,
#' # but does add potentially useful hclust dendrogram to the sides...
#' heatmap(otu_table(gpac))
#' plot_heatmap(gpac)
#' # example relabelling based on a sample variable and taxonomic rank.
#' plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family")
#' # Now repeat the plot, but change the color scheme in various ways.
#' # See the online tutorial for many other examples.
plot_heatmap <- function(physeq, method="NMDS", distance="bray", 
	sample.label=NULL, taxa.label=NULL, 
	low="#000033", high="#66CCFF", na.value="black", trans=log_trans(4), 
	max.label=250, title=NULL, species.label=taxa.label, ...){
	
	if( is.null(taxa.label) & !is.null(species.label) ){
		# If species.label was provided, but not taxa.label, assign it to taxa.label
		taxa.label = species.label
	}
	
	# Initialize sample and species order vectors as NULL
	OTUorder <- sample.order <- NULL
	
	# Copy the approach from NeatMap by doing ordination on samples, but use 
	# phyloseq-wrapped distance/ordination procedures.
	# Reorder by the angle in radial coordinates on the 2-axis plane.
	if( !is.null(method) ){
		# Capture the NMDS iterations cat() output with capture.output
		#junk = capture.output( ps.ord <- ordinate(physeq, method, distance), file=NULL)
		junk = capture.output( ps.ord <- ordinate(physeq, method, distance, ...), file=NULL)
		reduction.result = scores(ps.ord, choices=c(1, 2), display="sites")
		sample.order = sample_names(physeq)[order(RadialTheta(reduction.result))]

		# Also want to re-order species, if possible
		test <- try(scores(ps.ord, choices=c(1, 2), display="species"), TRUE)
		if( class(test) != "try-error" & !is.null(test) ){			
			OTUreduct = scores(ps.ord, choices=c(1, 2), display="species")
			OTUorder  = taxa_names(physeq)[order(RadialTheta(OTUreduct))]
		}
	}

	# melt physeq with the standard user-accessible data melting function
	# for creating plot-ready data.frames, psmelt.
	# This is also used inside some of the other plot_* functions.
	adf = psmelt(physeq)	
	# Coerce the main axis variables to character. 
	# Will reset them to factor if re-ordering is needed.
	adf$OTU = as(adf$OTU, "character")
	adf$Sample = as(adf$Sample, "character")
	if( !is.null(sample.order) ){
		# If sample-order is available, coerce to factor with special level-order
		adf$Sample = factor(adf$Sample, levels=sample.order)
	} else {
		# Make sure it is a factor, but with default order/levels
		adf$Sample = factor(adf$Sample)
	}
	if( !is.null(OTUorder) ){
		# If OTU-order is available, coerce to factor with special level-order
		adf$OTU = factor(adf$OTU, levels=OTUorder)
	} else {
		# Make sure it is a factor, but with default order/levels
		adf$OTU = factor(adf$OTU)
	}

	## Now the plotting part
	# Initialize p.
	p = ggplot(adf, aes(Sample, OTU, fill=Abundance)) + geom_tile()

	# # Don't render labels if more than max.label
	# Samples
	if( nsamples(physeq) <= max.label ){
		# Add resize layer for samples if there are fewer than max.label
		p <- p + theme(
			axis.text.x = element_text(
				size=manytextsize(nsamples(physeq), 4, 30, 12),
				angle=-90, vjust=0.5, hjust=0
			)
		)		
	} else {
		p = p + scale_x_discrete("Sample", labels="")
	}
	# OTUs
	if( ntaxa(physeq) <= max.label ){
		# Add resize layer for OTUs if there are fewer than max.label
		p <- p + theme(
			axis.text.y = element_text(
				size=manytextsize(ntaxa(physeq), 4, 30, 12)
			)
		)
	} else {
		# Remove the labels from any rendering.
		p = p + scale_y_discrete("OTU", labels="")
	}
	
	# # Axis Relabeling (Skipped if more than max.label):
	# Re-write sample-labels to some sample variable...
	if( !is.null(sample.label) & nsamples(physeq) <= max.label){
		# Make a sample-named char-vector of the values for sample.label
		labvec = as(get_variable(physeq, sample.label), "character")
		names(labvec) <- sample_names(physeq)
		if( !is.null(sample.order) ){
			# Re-order according to sample.order
			labvec = labvec[sample.order]			
		}
		# Replace any NA (missing) values with "" instead. Avoid recycling labels.
		labvec[is.na(labvec)] <- ""
		# Add the sample.label re-labeling layer
		p = p + scale_x_discrete(sample.label, labels=labvec)
	}
	if( !is.null(taxa.label) & ntaxa(physeq) <= max.label){
		# Make a OTU-named vector of the values for taxa.label
		labvec <- as(tax_table(physeq)[, taxa.label], "character")
		names(labvec) <- taxa_names(physeq)
		if( !is.null(OTUorder) ){		
			# Re-order according to OTUorder
			labvec <- labvec[OTUorder]
		}
		# Replace any NA (missing) values with "" instead. Avoid recycling labels.
		labvec[is.na(labvec)] <- ""		
		# Add the taxa.label re-labeling layer
		p = p + scale_y_discrete(taxa.label, labels=labvec)
	}
	
	# Color scale transformations
	if( !is.null(trans) ){
		p = p + scale_fill_gradient(low=low, high=high, trans=trans, na.value=na.value)
	} else {
		p = p + scale_fill_gradient(low=low, high=high, na.value=na.value)	
	}
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p = p + ggtitle(title)
	}
			
	return(p)
}
################################################################################
################################################################################