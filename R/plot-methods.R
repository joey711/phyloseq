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
		ape::plot.phylo(tree, ...)	
		ape::nodelabels(as.character(1:max(tree$edge)), node=1:max(tree$edge))
		ape::edgelabels(as.character(1:nrow(tree$edge)), edge=1:nrow(tree$edge))		
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
#' @importFrom reshape melt
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph get.edgelist
#' @export
#' @examples 
#' 
#' data(enterotype)
#' ig <- make_sample_network(enterotype, FALSE, max.dist=0.3)
#' plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # Change distance parameter
#' ig <- make_sample_network(enterotype, FALSE, max.dist=0.2)
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
	p <- p + ggplot2::theme_bw() + 
			ggplot2::opts(
				panel.grid.major = ggplot2::theme_blank(), 
				panel.grid.minor = ggplot2::theme_blank(), 
				axis.text.x      = ggplot2::theme_blank(),
				axis.text.y      = ggplot2::theme_blank(),
				axis.title.x     = ggplot2::theme_blank(),
				axis.title.y     = ggplot2::theme_blank(),
				axis.ticks       = ggplot2::theme_blank(),
				panel.border     = ggplot2::theme_blank()
			)

	# Add the graph vertices as points
	p <- p + ggplot2::geom_point(aes_string(color=color, shape=shape), size=point_size)

	# Add the text labels
	if( !is.null(label) ){
		p <- p + ggplot2::geom_text(ggplot2::aes_string(label=label), size = 2, hjust=hjust)		
	}
	
	# Add the edges:
	p <- p + ggplot2::geom_line(ggplot2::aes_string(group="id", color=line_color), 
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
#' @importFrom reshape melt
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
	richness_map <- ggplot2::aes_string(x=x, y="value", color=color, shape=shape)		
	
	# Make the ggplot. Note that because ggplot2 is fully loaded in namespace,
	# its functions must be fully-qualified (ggplot2::), according to Bioconductor rules
	p <- ggplot2::ggplot(mdf, richness_map) + 
		ggplot2::geom_point(size=2) + 
		ggplot2::geom_errorbar(ggplot2::aes(ymax=value + se, ymin=value - se), width=0.2) +	
		ggplot2::opts(axis.text.x = ggplot2::theme_text(angle = -90, hjust = 0)) +
		ggplot2::scale_y_continuous('richness [number of species]') +
		ggplot2::facet_grid(~variable) 
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
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{sampleData-class}}. The data about which you want to 
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
#' # # DPCoA - because scores.dpcoa() is internal to phyloseq, can only test
#' # # this chunk when new plot_ordination() is installed
#' # topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
#' # GP1   <- prune_species(topsp, GP)
#' # GP.dpcoa <- DPCoA(GP1)
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
		ord_map <- ggplot2::aes_string(x=x, y=y, color=color, shape=shape, na.rm=TRUE)
	} else if(type=="biplot"){
		# biplot, id.type must map to color or size. Only color if none specified.
		if( is.null(color) ){
			ord_map <- ggplot2::aes_string(x=x, y=y, color="id.type",
							shape=shape, na.rm=TRUE)
		} else {
			ord_map <- ggplot2::aes_string(x=x, y=y, size="id.type",
							color=color, shape=shape, na.rm=TRUE)
		}
	}

	# Plot-building section
	p <- ggplot2::ggplot(DF, ord_map) + ggplot2::geom_point(na.rm=TRUE)
	
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
		label_map <- ggplot2::aes_string(x=x, y=y, label=label, na.rm=TRUE)
		p <- p + ggplot2::geom_text(label_map, data=rm.na.phyloseq(DF, label),
					size=2, vjust=1.5, na.rm=TRUE)
	}

	if( !is.null(title) ){
		p <- p + ggplot2::opts(title = title)
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
#' # gpca <- cca.phyloseq(GP)
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