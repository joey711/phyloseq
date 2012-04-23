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
#'  \code{\link{plot_ordination}}
#'  \code{\link{plot_taxa_bar}}
#'  \code{\link{plot_sample_network}}
#'  \code{\link{plot_tree_phyloseq}}
#'  \code{\link{plot_richness_estimates}}
#'
#' @export
#' @docType methods
#' @rdname plot_phyloseq-methods
#'
#' @examples 
#'  ## data(esophagus)
#'  ## plot_phyloseq(esophagus)
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
		plot.phylo(tree, ...)	
		nodelabels(as.character(1:max(tree$edge)), node=1:max(tree$edge))
		edgelabels(as.character(1:nrow(tree$edge)), edge=1:nrow(tree$edge))		
	} else {
		plot_richness_estimates(physeq)
	}
})
################################################################################
################################################################################
#' Plot sample-wise microbiome network (ggplot2)
#'
#' A custom plotting function for displaying graph objects created by 
#' \code{\link[igraph]{igraph}} from a 
#' phylogenetic sequencing experiment (\code{\link{phyloseq-class}}),
#' using advanced \code{\link[ggplot2]{ggplot}}2 formatting.
#'
#' @usage plot_sample_network(g, physeq=NULL,
#' 	color=NULL, shape=NULL, point_size=4, alpha=1,
#' 	label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold)
#'
#' @param g (Required). An \code{\link[igraph]{igraph}}-class object created
#'  either by the convenience wrapper \code{\link{make_sample_network}}, 
#'  or directly by the tools in the igraph-package.
#'
#' @param physeq (Optional). Default \code{NULL}. 
#'  A \code{\link{phyloseq-class}} object on which \code{g} is based.
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
#'  for drawing a graph. Should be able to take an \code{\link{igraph}}-class
#'  as sole argument, and return a two-column coordinate matrix with \code{nrow}
#'  equal to the number of vertices. For possible options already included in 
#'  \code{igraph}-package, see the others also described in the help file:
#' 
#' \code{\link[igraph]{layout.fruchterman.reingold}}
#'
#' @return A \code{\link{ggplot}}2 plot.
#' 
#' @seealso 
#'  \code{\link{make_sample_network}}
#'
#' @references
#'  Code modified from code now hosted on GitHub by Scott Chamberlain:
#'  \url{https://github.com/SChamberlain/gggraph}
#'
#'  The code most directly used/modified was first posted here:
#'  \url{http://www.r-bloggers.com/basic-ggplot2-network-graphs/}
#' 
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom reshape melt.data.frame
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph get.edgelist
#' @export
#' @examples 
#' 
#' data(enterotype)
#' ig <- make_sample_network(enterotype, max.dist=0.3)
#' plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # Change distance parameter
#' ig <- make_sample_network(enterotype, max.dist=0.2)
#' plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
plot_sample_network <- function(g, physeq=NULL,
	color=NULL, shape=NULL, point_size=4, alpha=1,
	label="value", hjust = 1.35, 
	line_weight=0.5, line_color=color, line_alpha=0.4,
	layout.method=layout.fruchterman.reingold){

	# Make the edge-coordinates data.frame
	edgeDF    <- data.frame(get.edgelist(g))
	edgeDF$id <- 1:length(edgeDF[, 1])

	# Make the vertices-coordinates data.frame
	vertDF    <- layout.method(g)
	colnames(vertDF) <- c("x", "y")
	vertDF    <- data.frame(value=g[[9]][[3]][["name"]], vertDF)
	
	# If phyloseq object provided, add its sample data to vertDF
	if( !is.null(physeq) ){
		SD     <- sampleData(physeq)[as.character(vertDF$value), ]
		vertDF <- data.frame(vertDF, SD) 
	}

	# Combine vertex and edge coordinate data.frames
	graphDF   <- merge(melt(edgeDF, id="id"), vertDF, by = "value") 
 
	# Initialize the ggplot
	p <- ggplot(vertDF, aes(x, y)) 

	# Strip all the typical annotations from the plot, leave the legend
	p <- p + theme_bw() + 
			opts(
				panel.grid.major = theme_blank(), 
				panel.grid.minor = theme_blank(), 
				axis.text.x      = theme_blank(),
				axis.text.y      = theme_blank(),
				axis.title.x     = theme_blank(),
				axis.title.y     = theme_blank(),
				axis.ticks       = theme_blank(),
				panel.border     = theme_blank()
			)

	# Add the graph vertices as points
	p <- p + geom_point(aes_string(color=color, shape=shape), size=point_size)

	# Add the text labels
	if( !is.null(label) ){
		p <- p + geom_text(aes_string(label=label), size = 2, hjust=hjust)		
	}
	
	# Add the edges:
	p <- p + geom_line(aes_string(group="id", color=line_color), 
				graphDF, size=line_weight, alpha=line_alpha)
	
	return(p)
}
################################################################################
################################################################################
#' Plot richness estimates, flexibly with ggplot2
#'
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
#' @usage plot_richness_estimates(physeq, x, color=NULL, shape=NULL)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otuTable-class}}. The data about which you want to estimate
#'  the richness.
#'
#' @param x (Optional). A variable to map to the horizontal axis. The vertical
#'  axis will be mapped to richness estimates and have units of total species.
#'  This parameter (\code{x}) can be either a character string indicating a
#'  variable in \code{sampleData} 
#'  (among the set returned by \code{sample.variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'
#'  The default value is \code{"sample.names"}, which will map each sample's name
#'  to a separate horizontal position in the plot.
#'
#' @param color (Optional). Default \code{NULL}. The sample variable to map
#'  to different colors. Like \code{x}, this can be a single character string 
#'  of the variable name in 
#'  \code{sampleData} 
#'  (among the set returned by \code{sample.variables(physeq)} );
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
#'  \code{sampleData} 
#'  (among the set returned by \code{sample.variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'  The shape scale is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_shape_manual}}.
#'
#' @return A \code{\link{ggplot}} plot object summarizing
#'  the richness estimates, and their standard error.
#' 
#' @seealso 
#'  \code{\link{estimate_richness}},
#'  \code{\link[vegan]{estimateR}},
#'  \code{\link[vegan]{diversity}}
#'
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom reshape melt.data.frame
#' @export
#' @examples 
#' # data(GlobalPatterns)
#' # plot_richness_estimates(GlobalPatterns, "SampleType")
#' # plot_richness_estimates(GlobalPatterns, "SampleType", "SampleType")
#' #
#' # # Define a human-associated versus non-human categorical variable:
#' # GP <- GlobalPatterns
#' # human.levels <- levels( getVariable(GP, "SampleType") ) %in% 
#' # c("Feces", "Mock", "Skin", "Tongue")
#' # human <- human.levels[getVariable(GP, "SampleType")]
#' # names(human) <- sample.names(GP)
#' # # Replace current SD with new one that includes human variable:
#' # sampleData(GP) <- sampleData(data.frame(sampleData(GP), human))
#' # 
#' # # Can use new "human" variable within GP as a discrete variable in the plot
#' # plot_richness_estimates(GP, "human", "SampleType")
#' # plot_richness_estimates(GP, "SampleType", "human")
#' #
#' # # Can also provide custom factor directly:
#' # plot_richness_estimates(GP, "SampleType", human)
#' # plot_richness_estimates(GP, human, "SampleType")
#' # 
#' # # Not run: Should cause an error:
#' # plot_richness_estimates(GP, "value", "value")
#' # #
plot_richness_estimates <- function(physeq, x="sample.names", color=NULL, shape=NULL){	
	# Make the plotting data.frame 
	DF <- data.frame(estimate_richness(physeq), sampleData(physeq))
	
	# If there is no "sample.names" variable in DF, add it
	if( !"sample.names" %in% names(DF) ){
		DF <- data.frame(DF, sample.names=sample.names(physeq))		
	}

	# If manually-supplied x, color, shape, add to DF, dummy var_name
	if(length(x) > 1){
		DF$x <- x
		names(DF)[names(DF)=="x"] <- deparse(substitute(x))
		x <- deparse(substitute(x))
	}
	if(length(color) > 1){
		DF$color <- color
		names(DF)[names(DF)=="color"] <- deparse(substitute(color))
		color <- deparse(substitute(color))
	}
	if(length(shape) > 1){
		DF$shape <- shape
		names(DF)[names(DF)=="shape"] <- deparse(substitute(shape))
		shape <- deparse(substitute(shape))
	}
	
	# melt, for different estimates
	if( is.null(color) | identical(x, color) ){
		mdf <- melt(DF[, c("S.obs", "S.chao1", "S.ACE", x)], 
			id=c(x))
	} else {
		mdf <- melt(DF[, c("S.obs", "S.chao1", "S.ACE", x, color)], 
			id=c(x, color))			
	}
			
	# Add standard error to melted df
	mdf    <- data.frame(mdf, se = c(rep(NA, nrow(DF)), DF[, "se.chao1"], DF[, "se.ACE"]) )	
	
	# map variables
	richness_map <- aes_string(x=x, y="value", color=color, shape=shape)		
	
	# Make the ggplot.
	p <- ggplot(mdf, richness_map) + 
		geom_point(size=2) + 
		geom_errorbar(aes(ymax=value + se, ymin=value - se), width=0.2) +	
		opts(axis.text.x = theme_text(angle = -90, hjust = 0)) +
		scale_y_continuous('richness [number of species]') +
		facet_grid(~variable) 
	return(p)
}
################################################################################
################################################################################
# The general case, could plot samples, taxa, or both (biplot/split). Default samples.
################################################################################
#' General ordination plotter based on ggplot2.
#'
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
#'  of ordination are defined by \code{R} packages. The supported classes 
#'  should be listed explicitly, but in the meantime, all ordination classes
#'  currently supported by the \code{\link[vegan]{scores}} function are
#'  supported here. There is no default, as the expectation is that the 
#'  ordination will be performed and saved prior to calling this plot function.
#'
#' @param type (Optional). The plot type. Default is \code{"samples"}. The
#'  currently supported options are 
#'  \code{c("samples", "sites", "species", "taxa", "biplot", "split")}.
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
#'  (among the set returned by \code{sample.variables(physeq)} )
#'  or
#'  taxonomic rank
#'  (among the set returned by \code{rank.names(physeq)}).
#'  
#'  Alternatively, if \code{type} indicates a single-plot 
#'  (\code{"samples"} or \code{"species"}), then
#'  it is also possible to supply a custom vector with length equal to
#'  the relevant number of samples or species
#'  (\code{nsamples(physeq)} or \code{nspecies(physeq)}).
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
#' @param title (Optional). Default \code{NULL}. Character string. The
#'  title to include over the plot. 
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
#'  \code{\link{plot_phyloseq}}
#'
#' @import ggplot2
#' @export
#' @examples 
#' ##
#' # data(GlobalPatterns)
#' # # Define a human-associated versus non-human binary variable:
#' # human.levels <- levels( getVariable(GlobalPatterns, "SampleType") ) %in%
#' 		# c("Feces", "Mock", "Skin", "Tongue")
#' # human <- human.levels[getVariable(GlobalPatterns, "SampleType")]
#' # names(human) <- sample.names(GlobalPatterns)
#' # # Need to clean the zeros from GlobalPatterns:
#' # GP <- prune_species(speciesSums(GlobalPatterns)>0, GlobalPatterns)
#' # # Get the names of the most-abundant
#' # top.TaxaGroup <- sort(
#' 		# tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE),
#' 		# decreasing = TRUE)
#' # top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
#' # # Now prune further, to just the most-abundant phyla
#' # GP <- subset_species(GP, Phylum %in% names(top.TaxaGroup))
#' # topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
#' # GP1   <- prune_species(topsp, GP)
#' # GP.dpcoa <- ordinate(GP1, "DPCoA")
#' # plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum")
#' # plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType") + geom_line() + geom_point(size=5)
#' # plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType", shape=human) + 
#'      # geom_line() + geom_point(size=5)
#' # plot_ordination(GP1, GP.dpcoa, type="species", color="Phylum") + geom_line() + geom_point(size=5)
#' # plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum", label="SampleType")
#' # plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum")
#' # plot_ordination(GP1, GP.dpcoa, type="biplot", color="Phylum")
#' # plot_ordination(GP1, GP.dpcoa, type="biplot", label="Phylum")
#' # plot_ordination(GP1, GP.dpcoa, type="split", color="Phylum", label="SampleType")
#' # plot_ordination(GP1, GP.dpcoa, type="split", color="SampleType", shape="Phylum", label="SampleType")
plot_ordination <- function(physeq, ordination, type="samples", axes=c(1, 2),
	color=NULL, shape=NULL, label=NULL, title=NULL, justDF=FALSE){

	if(class(physeq)!="phyloseq"){stop("physeq must be phyloseq-class.")}
	if(type == "samples"){type <- "sites"} # Compatibility with phyloseq
	if(type == "taxa"){type <- "species"} # Compatibility with phyloseq
	if( !type %in% c("sites", "species", "biplot", "split") ){stop("type argument not supported.")}

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
		specDF$id.type <- "species"
		siteDF$id.type <- "samples"
		# Merge the two data.frame together, for joint plotting.
		DF <- merge(specDF, siteDF, all=TRUE)
		# Replace NA with "sample" or "species", where appropriate (factor/character)
		if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "samples") }
		if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "species") }
		if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "samples") }
		if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "species") }		
	}
	# In case user wants the plot-DF for some other purpose, return early
	if(justDF){return(DF)}
	
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
			p <- p + scale_color_discrete(name="type")
		} else {
			# Adjust size so that samples are bigger than species by default.
			p <- p + scale_size_manual("type", values=c(samples=5, species=2))		
		}
	}

	# Add the text labels
	if( !is.null(label) ){
		label_map <- aes_string(x=x, y=y, label=label, na.rm=TRUE)
		p <- p + geom_text(label_map, data=rm.na.phyloseq(DF, label),
					size=2, vjust=1.5, na.rm=TRUE)
	}

	if( !is.null(title) ){
		p <- p + opts(title = title)
	}
	return(p)
}
################################################################################
# Define the ord.plot.DF.internal
################################################################################
#' @keywords internal
ord.plot.DF.internal <- function(physeq, ordination, type="samples", axes=c(1, 2)){

	coord <- scores(ordination, choices=axes, display=type)
	# coord row.names index order should match physeq. Enforce.
	if( type == "species" ){
		coord <- coord[species.names(physeq), ]
	} else if(type == "sites"){
		coord <- coord[sample.names(physeq), ]		
	}
	
	# If there is supplemental data, add it, else, return coord
	supp <- NULL
	# Define supplemental data
	if( !is.null(sampleData(physeq, FALSE)) & type == "sites"){
		supp  <- sampleData(physeq) # Supplemental data, samples
	}else if( !is.null(taxTab(physeq, FALSE)) & type == "species"){
		supp  <- taxTab(physeq) # Supplemental data, taxa
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
#'  \code{\link{plot_ordination}}
#'
#' @import ggplot2
#' @export
#' @examples 
#' ##
#' # data(GlobalPatterns)
#' # # Need to clean the zeros from GlobalPatterns:
#' # GP <- prune_species(speciesSums(GlobalPatterns)>0, GlobalPatterns)
#' # sampleData(GP)$human <- factor(human)
#' # # Get the names of the most-abundant phyla
#' # top.TaxaGroup <- sort(
#' #   	tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE),
#' #   decreasing = TRUE)
#' # top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
#' # # Prune to just the most-abundant phyla
#' # GP <- subset_species(GP, Phylum %in% names(top.TaxaGroup))
#' # # Perform a correspondence analysis
#' # gpca <- ordinate(GP, "CCA")
#' # # # Make species topo with a subset of points layered
#' # # First, make a basic plot of just the species
#' # p1 <- plot_ordination(GP, gpca, "species", color="Phylum")
#' # # Re-draw this as topo without points, and facet
#' # p1 <- ggplot(p1$data, p1$mapping) + geom_density2d() + facet_wrap(~Phylum)
#' # # Add a layer of a subset of species-points that are furthest from origin.
#' # p53 <- p1 + geom_point(data=subset_ord_plot(p1, 1.0, "square"), size=1) 
#' # print(p53)
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
#' @usage plot_taxa_bar(otu, taxavec="Domain",
#'	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category,  
#'	facet_formula = . ~ TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE)
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
#' ##
#' # data(enterotype)
#' # TopNOTUs <- names(sort(speciesSums(enterotype), TRUE)[1:10]) 
#' # ent10   <- prune_species(TopNOTUs, enterotype)
#' # (p <- plot_taxa_bar(ent10, "Genus", x="SeqTech", fill="TaxaGroup") +
#' #    facet_wrap(~Enterotype) )
plot_taxa_bar <- function(otu, taxavec="Domain",
	showOnlyTheseTaxa=NULL, threshold=NULL, x_category="sample", fill_category=x_category, 
	facet_formula = . ~ TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE){

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

	########################################
	# Build the ggplot
	p  <- ggplot(df) + 
		opts(axis.text.x=theme_text(angle=-90, hjust=0))

	p <- p + 
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
		# Some reasonable default options
		opts(panel.grid.minor = theme_blank()) + 
		opts(panel.grid.major = theme_blank()) +
		opts(panel.border = theme_blank()) +
		labs(y="Relative Abundance", x=x_category, fill=fill_category)
		
	# Should the individual OTU points be added. Default FALSE
	if( OTUpoints ){
		p <- p + geom_point(
			data=df, 
			eval(call("aes",
				x=as.name(x_category), 
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
				x=as.name(x_category), 
				y=quote(Abundance+0.01), 
				label=quote(ID)
			)),
		)
	}
		

	if( !is.null(facet_formula) ){	
		p <- p + facet_grid(facet_formula)
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
#' ## data(GlobalPatterns)
#' ## # for reproducibility
#' ## set.seed(711)
#' ## ex2 <- prune_species(sample(species.names(GlobalPatterns), 50), GlobalPatterns)
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
#' Plot and annotate tree-tip using base/ape graphics
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
#' @import ape
#' @export
#'
#' @examples
#' # data(GlobalPatterns)
#' # GP <- GlobalPatterns
#' # GP.chl <- subset_species(GP, Phylum=="Chlamydiae")
#' # plot_tree_phyloseq(GP.chl, color_factor="SampleType",
#' # 			type_abundance_value=TRUE, 
#' # 			treeTitle="Chlamydiae in Global Patterns Data")
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
	plot.phylo(tree, type="phylogram", show.tip.label=FALSE,
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
#'
#' @keywords internal
tree.layout <- function(
  phylo,
  layout = 'default',
  layout.ancestors = FALSE,
  align.seq.names = NA
) {
  # Number of nodes and leaves.
  n.nodes <- length(phylo$tip.label)+phylo$Nnode
  n.leaves <- length(phylo$tip.label)

  t.labels <- phylo$tip.label
  n.labels <- ((n.leaves+1):n.nodes)
  if (!is.null(phylo$node.label)) {
    n.labels <- phylo$node.label
  }

  # Create the skeleton data frame.
  df <- data.frame(
                   node=c(1:n.nodes),                                            # Nodes with IDs 1 to N.
                   angle=0,
                   x=0,                                                          # These will contain the x and y coordinates after the layout procedure below.
                   y=0,
                   label=c(t.labels, n.labels),            # The first n.leaves nodes are the labeled tips.
                   is.leaf=c(rep(TRUE, n.leaves), rep(FALSE, n.nodes-n.leaves)),    # Just for convenience, store a boolean whether it's a leaf or not.
                   parent=0,                                                     # Will contain the ID of the current node's parent
                   children=0,                                                   # Will contain a list of IDs of the current node's children
                   branch.length=0                                               # Will contain the branch lengths
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
  df[df$is.leaf==TRUE,]$y = c(1:n.leaves)

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
treetextsize <- function(n){
	# empirically chosen size-value calculator.
	s <- 6 * exp(-n/100)
	# enforce a floor.
	s <- ifelse(s > 0.5, s, 0.5)
	# enforce a max
	s <- ifelse(s < 4, s, 4)
	return(s)
}
################################################################################
# Define an internal function for mapping phyloseq data variables to melted.tip
#' @keywords internal
treeMapVar2Tips <- function(melted.tip, physeq, variate){
	# If variate is taxTab-variable: Map taxTab-variable to melted.tip
	if( variate %in% rank.names(physeq, FALSE) ){
		# Add relevant taxTab column.
		x <- as(taxTab(physeq), "matrix")[, variate, drop=TRUE]
		names(x) <- species.names(physeq)
		return( x[as(melted.tip$species.names, "character")] )
	}
	# If variate is sampleMap-variable: Map sample-variable to melted.tip
	if( variate %in% sample.variables(physeq, FALSE) ){
		x <- as.vector(data.frame(sampleData(physeq))[, variate])
		names(x) <- sample.names(physeq)				
		return( x[as(melted.tip$variable, "character")] )
	}	
}
################################################################################
# The "tree only" setting. Simple. No annotations.
#' @keywords internal
plot_tree_only <- function(physeq){
	# Create the tree data.frame
	tdf <- tree.layout(tre(physeq))
	# build tree lines
	p <- ggplot(subset(tdf, type == "line")) + 
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) 
	return(p)
}
################################################################################
# The "sampledodge" plot_tree subset function.
#' @keywords internal
#' @importFrom scales log_trans
#' @importFrom reshape melt
#' @importFrom reshape melt.array
#' @importFrom plyr aaply
#' @importFrom plyr ddply
plot_tree_sampledodge <- function(physeq, color, shape, size, min.abundance, 
				label.tips, text.size, sizebase, base.spacing){

	# Create the tree data.frame
	tdf <- tree.layout(tre(physeq))
	
	# build tree lines
	p <- ggplot(subset(tdf, type == "line")) + 
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
								
	# Get the subset of tdf for just the tips (leaves)
	speciesDF <- subset(tdf, type=="label")
	
	# Add abundance data for each species
	# # First, re-order speciesDF ensure match with otuTable
	rownames(speciesDF) <- as(speciesDF$label, "character")
	speciesDF <- speciesDF[species.names(physeq), ]
	# # subset speciesDF to just what you need for tip plotting
	speciesDF <- data.frame(speciesDF[, c("x", "y")], species.names=rownames(speciesDF))
	
	# # Make the 0-values NA so they're not plotted. 
	OTU 		<- as(otuTable(physeq), "matrix") # Coerce to matrix.
	if(!speciesAreRows(physeq)){OTU <- t(OTU)} # Enforce orientation.
	OTU[OTU==0] <- NA
	
	# # Now add abundance table
	speciesDF 	<- data.frame(speciesDF, OTU)
	
	# # Now melt to just what you need for adding to plot
	melted.tip <- melt(speciesDF, id=c("x", "y", "species.names"))
	
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

	# Build the tip-label portion of the melted.tip data.frame, if needed.
	if( !is.null(label.tips) ){
		if( label.tips == "species.names" ){
			melted.tip$tipLabels <- melted.tip[, "species.names"]
		} else {
			melted.tip$tipLabels <- treeMapVar2Tips(melted.tip, physeq, label.tips)
		}
	} 	

	# color-map handling. Names "variable", "value" have specieal meaning.	
	if( !is.null(color) ){
		if( color %in% c("sample.names", "samples") ){
			color <- "variable"
		} else {
			melted.tip$color <- treeMapVar2Tips(melted.tip, physeq, color)
			names(melted.tip)[names(melted.tip)=="color"] <- color # rename to name of color variable
		}
	}

	# shape-map handling. Names "variable", "value" have specieal meaning.	
	if( !is.null(shape) ){
		if( shape %in% c("sample.names", "samples") ){
			shape <- "variable"
		} else if( !is.null(shape) ){
			melted.tip$shape <- treeMapVar2Tips(melted.tip, physeq, shape)
			names(melted.tip)[names(melted.tip)=="shape"] <- shape # rename to name of shape variable
		}
	}
	
	# size-map handling. Names "abundance", "variable", "value" have special meaning.
	if( !is.null(size) ){	
		if( size %in% c("abundance", "Abundance", "abund") ){
			size <- "value"
		} else if( !is.null(size) ){
			melted.tip$size <- treeMapVar2Tips(melted.tip, physeq, size)
			names(melted.tip)[names(melted.tip)=="size"] <- size # rename to name of size variable
		}
	}
		
	# The general tip-point map. Objects can be NULL, and that aesthetic gets ignored.
	tip.map <- aes_string(x="x + x.adj + x.spacer.base", y="y", color=color, shape=shape, size=size)
	
	# Add the new point layer.
	p <- p + geom_point(tip.map, data=melted.tip)

	# Optionally-add abundance value label to each point.
	# This size needs to match point size.
	if( any(melted.tip$value >= min.abundance[1]) ){
		if( is.null(size) ){
			point.label.map <- aes_string(x="x + x.adj + x.spacer.base", y="y", label="value")
			p <- p + geom_text( point.label.map, data=subset(melted.tip, value>=min.abundance[1]), size=1)
		} else {
			point.label.map <- aes_string(x="x + x.adj + x.spacer.base", y="y",
				label="value", size=paste("0.025*", size, sep=""))
			p <- p + geom_text( point.label.map, angle=45, hjust=0,
						data=subset(melted.tip, value>=min.abundance[1]) )
		}
	}

	# If no text.size given, calculate it from number of tips ("species", aka taxa)
	# This is very fast. No need to worry about whether text is printed or not. DRY.
	if( is.null(text.size) ){
		text.size <- treetextsize(nspecies(physeq))
	}

	# If indicated, add the species labels to the right of points.
	if( !is.null(label.tips) ){
		# melted.tip.far has only one row per tip,
		# the farthest horiz. adjusted position (one for each taxa)
		melted.tip.far <- ddply(melted.tip, "species.names", function(df){
			df[df$h.adj.index == max(df$h.adj.index), , drop=FALSE]
		})
		# Create the tip-label aesthetic map.
		label.map <- aes(x=x + x.adj + 2*x.spacer.base, y=y, label=tipLabels)
		# Add labels layer to plotting object.
		p <- p + geom_text(label.map, data=melted.tip.far, size=I(text.size), hjust=0)
	}
	
	# Adjust name of 
	if( !is.null(size) ){ 	
		# p <- p + scale_size_continuous("Abundance", trans=log_trans(sizebase))
		p <- p + scale_size_continuous(size, trans=log_trans(sizebase))
	}
	
	# Update legend-name of color or shape
	if( as.logical(sum(color == "variable")) ){
		p <- update_labels(p, list(colour = "Samples"))
	}
	if( as.logical(sum(shape == "variable")) ){
		p <- update_labels(p, list(shape  = "Samples"))
	}
			
	return(p)		
}
################################################################################
#' Plot a phylogenetic tree with optional annotations
#'
#' This function is intended to facilitate easy graphical investigation of 
#' the phylogenetic tree, as well as sample data. Note that for phylogenetic
#' sequencing of samples with large richness, some of the options in this 
#' function will be prohibitively slow to render, or too dense to be
#' interpretable. A rough ``rule of thumb'' is to use subsets of data 
#' with not many more than 200 taxa per plot, sometimes less depending on the
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
#' @usage plot_tree(physeq, method="sampledodge", color=NULL, shape=NULL, size=NULL,
#'  min.abundance=Inf, label.tips=NULL, text.size=NULL, sizebase=5, base.spacing = 0.02)
#'
#' @param physeq (Required). The data about which you want to 
#'  plot and annotate a phylogenetic tree, in the form of a
#'  single instance of the \code{\link{phyloseq-class}}, containing at 
#'  minimum a phylogenetic tree component (try \code{\link{tre}}).
#'  One of the major advantages of this function over basic tree-plotting utilities
#'  in the \code{\link{ape}}-package is the ability to easily annotate the tree
#'  with sample variables and taxonomic information. For these uses, 
#'  the \code{physeq} argument should also have a \code{\link{sampleData}}
#'  and/or \code{\link{taxTab}} component(s).
#' 
#' @param method (Optional). Character string. Default \code{"sampledodge"}. 
#'  The name of the annotation method to use. 
#'  This will be expanded in future versions.
#'  Currently only \code{"sampledodge"} and \code{"treeonly"} are supported.
#'  The \code{"sampledodge"} option results in points
#'  drawn next to leaves if individuals from that taxa were observed,
#'  and a separate point is drawn for each sample.
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
#'  If \code{"species.names"}, then the name of the taxa will be added 
#'  to the tree; either next to the leaves, or next to
#'  the set of points that label the leaves. Alternatively,
#'  if this is one of the rank names (from \code{rank.names(physeq)}),
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
#' @author Paul McMurdie, relying on supporting code from
#'  Gregory Jordan \email{gjuggler@@gmail.com}
#' 
#' @importFrom scales log_trans
#' @importFrom reshape melt
#' @export
#' @examples
#' # # # Using plot_tree with the esophagus dataset.
#' # data(esophagus)
#' # plot_tree(esophagus)
#' # plot_tree(esophagus, color="samples")
#' # plot_tree(esophagus, size="abundance")
#' # plot_tree(esophagus, size="abundance", color="samples")
#' # plot_tree(esophagus, size="abundance", color="samples", base.spacing=0.03)
#' # # #
#' # # # Using plot_tree with the Global Patterns dataset
#' # # Subset Global Patterns dataset to just the observed Archaea
#' # gpa <- subset_species(GlobalPatterns, Kingdom=="Archaea")
#' # # The number of different Archaeal species from this dataset is small enough ...
#' # nspecies(gpa)
#' # # ... that it is reasonable to consider displaying the phylogenetic tree directly.
#' # # (probably not true of the total dataset)
#' # nspecies(GlobalPatterns)
#' # # Some patterns are immediately discernable with minimal parameter choices:
#' # plot_tree(gpa, color="SampleType")
#' # plot_tree(gpa, color="Phylum")
#' # plot_tree(gpa, color="SampleType", shape="Phylum")
#' # plot_tree(gpa, color="Phylum", label.tips="Genus")
#' # # However, the text-label size scales with number of species, and with common
#' # # graphics-divice sizes/resolutions, these ~200 taxa still make for a 
#' # # somewhat crowded graphic. 
#' # # 
#' # # Let's instead subset further ot just the Crenarchaeota
#' # gpac <- subset_species(gpa, Phylum=="Crenarchaeota")
#' # plot_tree(gpac, color="SampleType", shape="Genus")
#' # plot_tree(gpac, color="SampleType", label.tips="Genus")
#' # # Let's add some abundance information.
#' # # Notice that the default spacing gets a little crowded when we map
#' # # species-abundance to point-size:
#' # plot_tree(gpac, color="SampleType", shape="Genus", size="abundance")
#' # # So let's spread it out a little bit with the base.spacing parameter.
#' # plot_tree(gpac, color="SampleType", shape="Genus", size="abundance", base.spacing=0.05)
plot_tree <- function(physeq, method="sampledodge", color=NULL, shape=NULL, size=NULL,
	min.abundance=Inf, label.tips=NULL, text.size=NULL,
	sizebase=5, base.spacing = 0.02){

	if( method %in% c("treeonly") ){
		p <- plot_tree_only(physeq)
	}
	
	if( method == "sampledodge" ){
		p <- plot_tree_sampledodge(physeq, color, shape, size, min.abundance, 
				label.tips, text.size, sizebase, base.spacing)
	}
	
	# Theme-ing:
	p <- p + opts(axis.ticks = theme_blank(),
			axis.title.x=theme_blank(), axis.text.x=theme_blank(),
			axis.title.y=theme_blank(), axis.text.y=theme_blank(),
			panel.background = theme_blank(),
			panel.grid.minor = theme_blank(),			
			panel.grid.major = theme_blank()
			)
	
	return(p)
}
################################################################################
################################################################################