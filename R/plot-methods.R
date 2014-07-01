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
#' \href{http://joey711.github.io/phyloseq}{phyloseq online tutorials}.
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
#'  \href{http://joey711.github.io/phyloseq/tutorials-index.html}{phyloseq frontpage tutorials}.
#'
#'  \code{\link{plot_ordination}}
#'  \code{\link{plot_heatmap}}
#'  \code{\link{plot_tree}}
#'  \code{\link{plot_network}}
#'  \code{\link{plot_bar}}
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
#' Microbiome Network Plot using ggplot2 
#'
#' There are many useful examples of phyloseq network graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_network-examples}{phyloseq online tutorials}.
#' A custom plotting function for displaying networks
#' using advanced \code{\link[ggplot2]{ggplot}}2 formatting.
#' The network itself should be represented using
#' the \code{igraph} package.
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
#' @param g (Required). An \code{igraph}-class object created
#'  either by the convenience wrapper \code{\link{make_network}}, 
#'  or directly by the tools in the igraph-package.
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
#'  for drawing a graph. Should be able to take an \code{igraph}-class
#'  as sole argument, and return a two-column coordinate matrix with \code{nrow}
#'  equal to the number of vertices. For possible options already included in 
#'  \code{igraph}-package, see the others also described in the help file:
#' 
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' \code{\link[igraph]{layout.fruchterman.reingold}}
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
#' @import reshape2
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph get.edgelist
#' @importFrom igraph get.vertex.attribute
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

  if( vcount(g) < 2 ){
    # Report a warning if the graph is empty
    stop("The graph you provided, `g`, has too few vertices. 
         Check your graph, or the output of `make_network` and try again.")
  }
  
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
	vertDF    <- data.frame(value=get.vertex.attribute(g, "name"), vertDF)
	
	# If phyloseq object provided,
	# AND it has the relevant additional data
	# THEN add it to vertDF
	if( !is.null(physeq) ){
		extraData <- NULL
		if( type == "samples" & !is.null(sample_data(physeq, FALSE)) ){
			extraData = data.frame(sample_data(physeq))[as.character(vertDF$value), , drop=FALSE]
		} else if( type == "taxa" & !is.null(tax_table(physeq, FALSE)) ){
			extraData =   data.frame(tax_table(physeq))[as.character(vertDF$value), , drop=FALSE]
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
#' Microbiome Network Plot using ggplot2 
#'
#' There are many useful examples of phyloseq network graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_net-examples}{phyloseq online tutorials}.
#' A custom plotting function for displaying networks
#' using advanced \code{\link[ggplot2]{ggplot}}2 formatting.
#' Note that this function is a performance and interface revision to
#' \code{\link{plot_network}}, which requires an \code{\link[igraph]{igraph}}
#' object as its first argument.
#' This new function is more in-line with other
#' \code{plot_*} functions in the \code{\link{phyloseq-package}}, in that its
#' first/main argument is a \code{\link{phyloseq-class}} instance.
#' Edges in the network are created if the distance between
#' nodes is below a (potentially arbitrary) threshold,
#' and special care should be given to considering the choice of this threshold.
#' However, network line thickness and opacity is scaled according to the
#' similarity of vertices (either samples or taxa),
#' helping to temper, somewhat, the effect of the threshold.
#' Also note that the choice of network layout algorithm can have a large effect
#' on the impression and interpretability of the network graphic,
#' and you may want to familiarize yourself with some of these options
#' (see the \code{laymeth} argument).
#'
#' @param physeq (Required). 
#'  The \code{\link{phyloseq-class}} object that you want to represent as a network.
#'  
#' @param distance (Optional). Default is \code{"bray"}. 
#'  Can be either a distance method supported by \code{\link[phyloseq]{distance}},
#'  or an already-computed \code{\link{dist}}-class with labels that match
#'  the indices implied by both the \code{physeq} and \code{type} arguments
#'  (that is, either sample or taxa names).
#'  If you used \code{\link[phyloseq]{distance}} to pre-calculate your \code{\link{dist}}ance,
#'  and the same \code{type} argument as provided here, then they will match.
#'  
#' @param maxdist (Optional). Default \code{0.7}. 
#'  The maximum distance value between two vertices
#'  to connect with an edge in the graphic.
#'
#' @param type (Optional). Default \code{"samples"}.
#'  Whether the network represented in the primary argument, \code{g},
#'  is samples or taxa/OTUs.
#'  Supported arguments are \code{"samples"}, \code{"taxa"},
#'  where \code{"taxa"} indicates using the taxa indices,
#'  whether they actually represent species or some other taxonomic rank.
#'  
#' @param laymeth (Optional). Default \code{"fruchterman.reingold"}.
#'  A character string that indicates the method that will determine
#'  the placement of vertices, typically based on conectedness of vertices
#'  and the number of vertices.
#'  This is an interesting topic, and there are lots of options.
#'  See \code{\link{igraph-package}} for related topics in general, 
#'  and see \code{\link[igraph]{layout.auto}} for descriptions of various
#'  alternative layout method options supported here.
#'  The character string argument should match exactly the
#'  layout function name with the \code{"layout."} omitted.
#'  Try \code{laymeth="list"} to see a list of options.
#'
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#' 
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).
#'  
#' @param rescale (Optional). Logical. Default \code{FALSE}.
#'  Whether to rescale the distance values to be \code{[0, 1]}, in which the
#'  min value is close to zero and the max value is 1.
#' 
#' @param point_size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#' 
#' @param point_alpha (Optional). Default \code{1}.
#'  A value between 0 and 1 for the alpha transparency of the vertex points.
#' 
#' @param point_label (Optional). Default \code{NULL}.
#'  The variable name in \code{physeq} covariate data to map to vertex labels.
#' 
#' @param hjust (Optional). Default \code{1.35}.
#'  The amount of horizontal justification to use for each label.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'  
#' @return A \code{\link{ggplot}}2 network plot.
#'  Will render to default graphic device automatically as print side effect.
#'  Can also be saved, further manipulated, or rendered to
#'  a vector or raster file using \code{\link{ggsave}}.
#' 
#' @seealso 
#'  Original network plotting functions:
#' 
#'  \code{\link{make_network}}
#' 
#'  \code{\link{plot_network}}
#' 
#' @import ggplot2
#' @import reshape2
#' @importFrom data.table data.table
#' @importFrom igraph layout.auto
#' @importFrom igraph layout.random
#' @importFrom igraph layout.circle
#' @importFrom igraph layout.sphere
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph layout.kamada.kawai
#' @importFrom igraph layout.spring
#' @importFrom igraph layout.reingold.tilford
#' @importFrom igraph layout.fruchterman.reingold.grid
#' @importFrom igraph layout.lgl
#' @importFrom igraph layout.graphopt
#' @importFrom igraph layout.svd
#' @importFrom igraph graph.data.frame
#' @importFrom igraph get.vertex.attribute
#' @export
#' @examples 
#' data(enterotype)
#' plot_net(enterotype, color="SeqTech", maxdist = 0.3)
#' plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "auto")
#' plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "svd")
#' plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "circle")
#' plot_net(enterotype, color="SeqTech", shape="Enterotype", maxdist = 0.3, laymeth = "circle")
plot_net <- function(physeq, distance="bray", type="samples", maxdist = 0.7,
                     laymeth="fruchterman.reingold", color=NULL, shape=NULL, rescale=FALSE,
                     point_size=5, point_alpha=1, point_label=NULL, hjust = 1.35, title=NULL){
  # Supported layout methods
  available_layouts = list(
    auto = layout.auto,
    random = layout.random,
    circle = layout.circle,
    sphere = layout.sphere,
    fruchterman.reingold = layout.fruchterman.reingold,
    kamada.kawai = layout.kamada.kawai,
    spring = layout.spring,
    reingold.tilford = layout.reingold.tilford,
    fruchterman.reingold.grid = layout.fruchterman.reingold.grid,
    lgl = layout.lgl,
    graphopt = layout.graphopt,
    svd = layout.svd
  )
  if(laymeth=="list"){
    return(names(available_layouts))
  }
  if(!laymeth %in% names(available_layouts)){
    stop("Unsupported argument to `laymeth` option. Please use an option returned by `plot_net(laymeth='list')`")
  }
  # 1. 
  # Calculate Distance
  if( inherits(distance, "dist") ){
    # If distance a distance object, use it rather than re-calculate
    Distance <- distance
    # Check that it at least has (a subset of) the correct labels
    possibleVertexLabels = switch(type, taxa=taxa_names(physeq), samples=sample_names(physeq))
    if( !all(attributes(distance)$Labels %in% possibleVertexLabels) ){
      stop("Some or all `distance` index labels do not match ", type, " names in `physeq`")
    }
  } else {
    # Coerce to character and attempt distance calculation
    scaled_distance = function(physeq, method, type, rescale=TRUE){
      Dist = distance(physeq, method, type)
      if(rescale){
        # rescale the distance matrix to be [0, 1]
        Dist <- Dist / max(Dist, na.rm=TRUE)
        Dist <- Dist - min(Dist, na.rm=TRUE)
      }
      return(Dist)
    }
    distance <- as(distance[1], "character")
    Distance = scaled_distance(physeq, distance, type, rescale)
  }
  # 2.
  # Create edge data.table
  dist_to_edge_table = function(Dist, MaxDistance=NULL, vnames = c("v1", "v2")){
    dmat <- as.matrix(Dist)
    # Set duplicate entries and self-links to Inf
    dmat[upper.tri(dmat, diag = TRUE)] <- Inf
    LinksData = data.table(reshape2::melt(dmat, varnames=vnames, as.is = TRUE))
    setnames(LinksData, old = "value", new = "Distance")
    # Remove self-links and duplicate links
    LinksData <- LinksData[is.finite(Distance), ]
    # Remove entries above the threshold, MaxDistance
    if(!is.null(MaxDistance)){
      LinksData <- LinksData[Distance < MaxDistance, ]
    }
    return(LinksData)
  }
  LinksData0 = dist_to_edge_table(Distance, maxdist)
  # 3. Create vertex layout
  # Make the vertices-coordinates data.table
  vertex_layout = function(LinksData, physeq=NULL, type="samples",
                           laymeth=igraph::layout.fruchterman.reingold, ...){
    # `physeq` can be anything, only has effect when non-NULL returned by sample_data or tax_table
    g = igraph::graph.data.frame(LinksData, directed=FALSE)
    vertexDT = data.table(laymeth(g, ...),
                          vertex=get.vertex.attribute(g, "name"))
    setkey(vertexDT, vertex)
    setnames(vertexDT, old = c(1, 2), new = c("x", "y"))
    extraData = NULL
    if( type == "samples" & !is.null(sample_data(physeq, FALSE)) ){
      extraData <- data.table(data.frame(sample_data(physeq)), key = "rn", keep.rownames = TRUE)
    } else if( type == "taxa" & !is.null(tax_table(physeq, FALSE)) ){
      extraData <- data.table(as(tax_table(physeq), "matrix"), key = "rn", keep.rownames = TRUE)
    }
    # Only mod vertexDT if extraData exists
    if(!is.null(extraData)){
      # Join vertexDT, extraData using data.table syntax. Presumes `vertex` is key in both.
      setnames(extraData, old = "rn", new = "vertex")
      vertexDT <- vertexDT[extraData]
      vertexDT <- vertexDT[!is.na(x), ]
    }
    return(vertexDT)
  }
  vertexDT = vertex_layout(LinksData0, physeq, type, available_layouts[[laymeth]])
  # 4.
  # Update the links layout for ggplot: x, y, xend, yend
  link_layout = function(LinksData, vertexDT){
    linkstart = vertexDT[LinksData$v1, x, y]
    linkend = vertexDT[LinksData$v2, x, y]
    setnames(linkend, old = c("y", "x"), new = c("yend", "xend"))
    LinksData <- cbind(LinksData, linkstart, linkend)
    return(LinksData)  
  }
  LinksData = link_layout(LinksData0, vertexDT)
  # 5.
  # Define ggplot2 network plot
  links_to_ggplot = function(LinksData, vertexDT, vertmap=aes(x, y)){
    p0 = ggplot(data=LinksData) + 
      geom_segment(aes(x, y, xend=xend, yend=yend, size=Distance, alpha=Distance)) +
      geom_point(mapping = vertmap, data=vertexDT, size=5, na.rm = TRUE) +
      scale_alpha(range = c(1, 0.1)) + 
      scale_size(range = c(2, 0.25))
    return(p0)
  }
  p = links_to_ggplot(LinksData, vertexDT,
                      vertmap = aes_string(x="x", y="y", color=color, shape=shape))
  # Add labels
  if(!is.null(point_label)){
    p <- p + geom_text(aes_string(x="x", y="y", label=point_label),
                       data = vertexDT, size = 2, hjust = hjust, na.rm = TRUE)
  }
  # Add default theme
  net_theme = theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x      = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.ticks       = element_blank(),
    panel.border     = element_blank()
  )
  p <- p + theme_bw() + net_theme
  return(p)
}
################################################################################
################################################################################
#' Plot alpha diversity, flexibly with ggplot2
#'
#' There are many useful examples of alpha-diversity graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_richness-examples}{phyloseq online tutorials}.
#' This function estimates a number of alpha-diversity metrics using the 
#' \code{\link{estimate_richness}} function,
#' and returns a \code{ggplot} plotting object. 
#' The plot generated by this function will include every sample
#' in \code{physeq}, but they can be further grouped on the horizontal axis
#' through the argument to \code{x}, 
#' and shaded according to the argument to \code{color} (see below).
#' You must use untrimmed, non-normalized count data for meaningful results, 
#' as many of these estimates are highly dependent on the number of singletons.
#' You can always trim the data later on if needed,
#' just not before using this function.
#'
#'  NOTE: Because this plotting function incorporates the output from 
#'  \code{\link{estimate_richness}}, the variable names of that output should
#'  not be used as \code{x} or \code{color} (even if it works, the resulting
#'  plot might be kindof strange, and not the intended behavior of this function).
#'  The following are the names you will want to avoid using in \code{x} or \code{color}:
#'
#'  \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otu_table-class}}. The data about which you want to estimate.
#'
#' @param x (Optional). A variable to map to the horizontal axis. The vertical
#'  axis will be mapped to the alpha diversity index/estimate
#'  and have units of total taxa, and/or index value (dimensionless).
#'  This parameter (\code{x}) can be either a character string indicating a
#'  variable in \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'
#'  The default value is \code{"samples"}, which will map each sample's name
#'  to a separate horizontal position in the plot.
#'
#' @param color (Optional). Default \code{NULL}. 
#'  The sample variable to map to different colors.
#'  Like \code{x}, this can be a single character string of the variable name in 
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
#' @param scales (Optional). Default \code{"free_y"}.
#'  Whether to let vertical axis have free scale that adjusts to
#'  the data in each panel.
#'  This argument is passed to \code{\link[ggplot2]{facet_wrap}}.
#'  If set to \code{"fixed"}, a single vertical scale will
#'  be used in all panels. This can obscure values if the
#'  \code{measures} argument includes both 
#'  richness estimates and diversity indices, for example.
#'  
#' @param nrow (Optional). Default is \code{1},
#'  meaning that all plot panels will be placed in a single row,
#'  side-by-side. 
#'  This argument is passed to \code{\link[ggplot2]{facet_wrap}}.
#'  If \code{NULL}, the number of rows and columns will be 
#'  chosen automatically (wrapped) based on the number of panels
#'  and the size of the graphics device.
#'
#' @param shsi (Deprecated). No longer supported. Instead see `measures` below.
#' 
#' @param measures (Optional). Default is \code{NULL}, meaning that
#'  all available alpha-diversity measures will be included in plot panels.
#'  Alternatively, you can specify one or more measures
#'  as a character vector of measure names.
#'  Values must be among those supported:
#'  \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#'
#' @param sortby (Optional). A character string subset of \code{measures} argument.
#'  Sort x-indices by the mean of one or more \code{measures},
#'  if x-axis is mapped to a discrete variable.
#'  Default is \code{NULL}, implying that a discrete-value horizontal axis
#'  will use default sorting, usually alphabetic. 
#'
#' @return A \code{\link{ggplot}} plot object summarizing
#'  the richness estimates, and their standard error.
#' 
#' @seealso 
#'  \code{\link{estimate_richness}}
#'  
#'  \code{\link[vegan]{estimateR}}
#'  
#'  \code{\link[vegan]{diversity}}
#'
#' There are many more interesting examples at the
#' \href{http://joey711.github.io/phyloseq/plot_richness-examples}{phyloseq online tutorials}.
#'
#' @import ggplot2
#' @import reshape2
#' @importFrom plyr is.discrete
#' @export
#' @examples 
#' ## There are many more interesting examples at the phyloseq online tutorials.
#' ## http://joey711.github.io/phyloseq/plot_richness-examples
#' data("soilrep")
#' plot_richness(soilrep, measures=c("InvSimpson", "Fisher"))
#' plot_richness(soilrep, "Treatment", "warmed", measures=c("Chao1", "ACE", "InvSimpson"), nrow=3)
#' data("GlobalPatterns")
#' plot_richness(GlobalPatterns, x="SampleType", measures=c("InvSimpson"))
#' plot_richness(GlobalPatterns, x="SampleType", measures=c("Chao1", "ACE", "InvSimpson"), nrow=3)
#' plot_richness(GlobalPatterns, x="SampleType", measures=c("Chao1", "ACE", "InvSimpson"), nrow=3, sortby = "Chao1")
plot_richness = function(physeq, x="samples", color=NULL, shape=NULL, title=NULL,
                          scales="free_y", nrow=1, shsi=NULL, measures=NULL, sortby=NULL){ 
  # Calculate the relevant alpha-diversity measures
  erDF = estimate_richness(physeq, split=TRUE, measures=measures)
  # Measures may have been renamed in `erDF`. Replace it with the name from erDF
  measures = colnames(erDF)
  # Define "measure" variables and s.e. labels, for melting.
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  # Remove any S.E. from `measures`
  measures = measures[!measures %in% ses]
	# Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(sample_data(physeq, errorIfNULL=FALSE)) ){
    # Include the sample data, if it is there.
	  DF <- data.frame(erDF, sample_data(physeq))
  } else {
    # If no sample data, leave it out.
    DF <- data.frame(erDF)
  }
	if( !"samples" %in% colnames(DF) ){
	  # If there is no "samples" variable in DF, add it
		DF$samples <- sample_names(physeq)
	}
	# sample_names used to be default, and should also work.
	# #backwardcompatibility
	if( !is.null(x) ){
		if( x %in% c("sample", "samples", "sample_names", "sample.names") ){
			x <- "samples"
		}
	} else {
    # If x was NULL for some reason, set it to "samples"
	  x <- "samples"
	}
	# melt to display different alpha-measures separately
	mdf = melt(DF, measure.vars=measures)
  # Initialize the se column. Helpful even if not used.
  mdf$se <- NA_integer_
  if( length(ses) > 0 ){
    ## Merge s.e. into one "se" column
    # Define conversion vector, `selabs`
    selabs = ses
    # Trim the "se." from the names
    names(selabs) <- substr(selabs, 4, 100)
    # Make first letter of selabs' names uppercase
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 1, 1))
    # use selabs conversion vector to process `mdf`
    mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs){selabs[i]}, selabs)
    for( i in 1:nrow(mdf) ){
      if( !is.na(mdf[i, "wse"]) ){
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    # prune the redundant columns
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  ## Interpret measures
  # If not provided (default), keep all 
  if( !is.null(measures) ){
    if( any(measures %in% as.character(mdf$variable)) ){
      # If any measures were in mdf, then subset to just those.
      mdf <- mdf[as.character(mdf$variable) %in% measures, ]
    } else {
      # Else, print warning about bad option choice for measures, keeping all.
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
	if( !is.null(shsi) ){
    # Deprecated:
    # If shsi is anything but NULL, print a warning about its being deprecated
		warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
	}
  # Address `sortby` argument
  if(!is.null(sortby)){
    if(!all(sortby %in% levels(mdf$variable))){
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if(!is.discrete(mdf[, x])){
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if(all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, x])){
      # Replace x-factor with same factor that has levels re-ordered according to `sortby`
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x],
                         levels = names(sort(tapply(X = mdf[wh.sortby, "value"],
                                                    INDEX = mdf[wh.sortby, x],
                                                    mean,
                                                    na.rm=TRUE, simplify = TRUE))))
    }
  }
  # Define variable mapping
  richness_map = aes_string(x=x, y="value", colour=color, shape=shape)
  # Make the ggplot.
  p = ggplot(mdf, richness_map) + geom_point(na.rm=TRUE)  
  # Add error bars if mdf$se is not all NA
  if( any(!is.na(mdf[, "se"])) ){
    p = p + geom_errorbar(aes(ymax=value + se, ymin=value - se), width=0.1) 
  }
  # Rotate horizontal axis labels, and adjust
	p = p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))
	# Add y-label 
	p = p + ylab('Alpha Diversity Measure') 
  # Facet wrap using user-options
	p = p + facet_wrap(~variable, nrow=nrow, scales=scales)
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
#' \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{phyloseq online tutorials}.
#' Convenience wrapper for plotting ordination results as a 
#' \code{ggplot2}-graphic, including
#' additional annotation in the form of shading, shape, and/or labels of
#' sample variables.
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
#'  Note that the color scheme is chosen automatically
#'  by \code{link{ggplot}},
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
#'  Many more examples are included in the 
#'  \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{phyloseq online tutorials}.
#'
#' Also see the general wrapping function:
#'
#' \code{\link{plot_phyloseq}}
#'
#' @import ggplot2
#' @importFrom vegan wascores
#' @export
#' @examples 
#' # See other examples at
#' # http://joey711.github.io/phyloseq/plot_ordination-examples
#' data(GlobalPatterns)
#' GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]), GlobalPatterns)
#' gp_bray_pcoa = ordinate(GP, "CCA", "bray")
#' plot_ordination(GP, gp_bray_pcoa, "samples", color="SampleType")
plot_ordination = function(physeq, ordination, type="samples", axes=1:2,
                            color=NULL, shape=NULL, label=NULL, title=NULL, justDF=FALSE){
  if(length(type) > 1){
    warning("`type` can only be a single option, but more than one provided. Using only the first.")
    type <- type[[1]]
  }
  if(length(color) > 1){
    warning("The `color` variable argument should have length equal to 1.",
            "Taking first value.")
    color = color[[1]][1]
  }
  if(length(shape) > 1){
    warning("The `shape` variable argument should have length equal to 1.",
            "Taking first value.")
    shape = shape[[1]][1]
  }
  if(length(label) > 1){
    warning("The `label` variable argument should have length equal to 1.",
            "Taking first value.")
    label = label[[1]][1]
  }
  official_types = c("sites", "species", "biplot", "split", "scree")
  if(!inherits(physeq, "phyloseq")){
    if(inherits(physeq, "character")){
      if(physeq=="list"){
        return(official_types)
      }
    } 
    warning("Full functionality requires `physeq` be phyloseq-class ",
            "with multiple components.")
  }
  # Catch typos and synonyms
  type = gsub("^.*site[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*species.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*taxa.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case=TRUE)
  type = gsub("^.*split[s]*.*$", "split", type, ignore.case=TRUE)
  type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case=TRUE)
  # If type argument is not supported...
  if( !type %in% official_types ){
    warning("type argument not supported. `type` set to 'samples'.\n",
            "See `plot_ordination('list')`")
    type <- "sites"
  }
  if( type %in% c("scree") ){
    # Stop early by passing to plot_scree() if "scree" was chosen as a type
    return( plot_scree(ordination, title=title) )
  }
  # Define a function to check if a data.frame is empty
  is_empty = function(x){
    length(x) < 2 | suppressWarnings(all(is.na(x)))
  }
  # The plotting data frames.
  # Call scores to get coordinates.
  # Silently returns only the coordinate systems available.
  # e.g. sites-only, even if species requested.
  specDF = siteDF = NULL
  try({siteDF <- scores(ordination, choices = axes, display="sites", physeq=physeq)}, silent = TRUE)
  try({specDF <- scores(ordination, choices = axes, display="species", physeq=physeq)}, silent = TRUE)
  # Check that have assigned coordinates to the correct object
  siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
  siteTaxaIntx = length(intersect(rownames(siteDF), taxa_names(physeq)))
  specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
  specTaxaIntx = length(intersect(rownames(specDF), taxa_names(physeq)))
  if(siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx){
    # Double-swap
    co = specDF
    specDF <- siteDF
    siteDF <- co
    rm(co)
  } else {
    if(siteSampIntx < specSampIntx){
      # Single swap
      siteDF <- specDF
      specDF <- NULL
    }
    if(specTaxaIntx < siteTaxaIntx){
      # Single swap 
      specDF <- siteDF
      siteDF <- NULL
    }
  }
  # If both empty, warn and return NULL
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  # If either is missing, do weighted average
  if(is_empty(specDF) & type != "sites"){
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)), stringsAsFactors=FALSE)
  }
  if(is_empty(siteDF) & type != "species"){ 
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))), stringsAsFactors=FALSE)
  }
  # Double-check that have assigned coordinates to the correct object
  specTaxaIntx <- siteSampIntx <- NULL
  siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
  specTaxaIntx <- length(intersect(rownames(specDF), taxa_names(physeq)))
  if(siteSampIntx < 1L & !is_empty(siteDF)){
    # If siteDF is not empty, but it doesn't intersect the sample_names in physeq, warn and set to NULL
    warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    siteDF <- NULL
  }
  if(specTaxaIntx < 1L & !is_empty(specDF)){
    # If specDF is not empty, but it doesn't intersect the taxa_names in physeq, warn and set to NULL
    warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    specDF <- NULL
  }
  # If you made it this far and both NULL, return NULL and throw a warning
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  if(type %in% c("biplot", "split") & (is_empty(siteDF) | is_empty(specDF)) ){
    # biplot and split require both coordinates systems available. 
    # Both were attempted, or even evaluated by weighted average.
    # If still empty, warn and switch to relevant type.
    if(is_empty(siteDF)){
      warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
      type <- "species"
    }
    if(is_empty(specDF)){
      warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
      type <- "sites"
    }
  }
  if(type != "species"){
    # samples covariate data frame, `sdf`
    sdf = NULL
    sdf = data.frame(access(physeq, slot="sam_data"), stringsAsFactors=FALSE)
    if( !is_empty(sdf) & !is_empty(siteDF) ){
      # The first two axes should always be x and y, the ordination axes.
      siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
    }
  }
  if(type != "sites"){
    # taxonomy data frame `tdf`
    tdf = NULL
    tdf = data.frame(access(physeq, slot="tax_table"), stringsAsFactors=FALSE)
    if( !is_empty(tdf) & !is_empty(specDF) ){
      # The first two axes should always be x and y, the ordination axes.
      specDF = cbind(specDF, tdf[rownames(specDF), ])
    }
  }
  # In "naked" OTU-table cases, `siteDF` or `specDF` could be matrix.
  if(!inherits(siteDF, "data.frame")){
    #warning("Sample Co-variables apparently missing in provided `physeq` for this plot-type. Coercing coord matrix to data.frame.")
    siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
  }  
  if(!inherits(specDF, "data.frame")){
    #warning("Taxonomy apparently missing in provided `physeq` for this plot-type. Coercing coord matrix to data.frame.")
    specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
  }
  # Define the main plot data frame, `DF`
  DF = NULL
  DF <- switch(EXPR = type, sites = siteDF, species = specDF, {
    # Anything else. In practice, type should be "biplot" or "split" here.
    # Add id.type label
    specDF$id.type <- "Taxa"
    siteDF$id.type <- "Samples"
    # But what if the axis variables differ b/w them?
    # Coerce specDF to match samples (siteDF) axis names
    colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
    # Merge the two data frames together for joint plotting.
    DF = merge(specDF, siteDF, all=TRUE)
    # Replace NA with "samples" or "taxa", where appropriate (factor/character)
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Samples") }
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Taxa") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Samples") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Taxa") }
    DF
  })
  # In case user wants the plot-DF for some other purpose, return early
  if(justDF){return(DF)}
  # Check variable availability before defining mapping.
  if(!is.null(color)){ 
    if(!color %in% names(DF)){
      warning("Color variable was not found in the available data you provided.",
              "No color mapped.")
      color <- NULL
    }
  }
  if(!is.null(shape)){ 
    if(!shape %in% names(DF)){
      warning("Shape variable was not found in the available data you provided.",
              "No shape mapped.")
      shape <- NULL
    }
  }
  if(!is.null(label)){ 
    if(!label %in% names(DF)){
      warning("Label variable was not found in the available data you provided.",
              "No label mapped.")
      label <- NULL
    }
  }
  # Grab the ordination axis names from the plot data frame (as strings)
  x = colnames(DF)[1]
  y = colnames(DF)[2]   
  # Mapping section
  if( ncol(DF) <= 2){
    # If there is nothing to map, enforce simple mapping.
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x=x, y=y)
  } else if( type %in% c("sites", "species", "split") ){
    ord_map = aes_string(x=x, y=y, color=color, shape=shape, na.rm=TRUE)
  } else if(type=="biplot"){
    # biplot, `id.type` should try to map to color and size. Only size if color specified.
    if( is.null(color) ){
      ord_map = aes_string(x=x, y=y, size="id.type", color="id.type", shape=shape, na.rm=TRUE)
    } else {
      ord_map = aes_string(x=x, y=y, size="id.type", color=color, shape=shape, na.rm=TRUE)
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
      p <- update_labels(p, list(colour="Ordination Type")) 
    } 
    # Adjust size so that samples are bigger than taxa by default.
    p <- p + scale_size_manual("type", values=c(Samples=5, Taxa=2))
  }
  # Add text labels to points
  if( !is.null(label) ){
    label_map <- aes_string(x=x, y=y, label=label, na.rm=TRUE)
    p = p + geom_text(label_map, data=rm.na.phyloseq(DF, label),
                      size=2, vjust=1.5, na.rm=TRUE)
  }
  # Optionally add a title to the plot
  if( !is.null(title) ){
    p = p + ggtitle(title)
  }
  # Add fraction variability to axis labels, if available
  if( length(extract_eigenvalue(ordination)[axes]) > 0 ){
    # Only attempt to add fraction variability
    # if extract_eigenvalue returns something
    eigvec = extract_eigenvalue(ordination)
    # Fraction variability, fracvar
    fracvar = eigvec[axes] / sum(eigvec)
    # Percent variability, percvar
    percvar = round(100*fracvar, 1)
    # The string to add to each axis label, strivar
    # Start with the curent axis labels in the plot
    strivar = as(c(p$label$x, p$label$y), "character")
    # paste the percent variability string at the end
    strivar = paste0(strivar, "   [", percvar, "%]")
    # Update the x-label and y-label
    p = p + xlab(strivar[1]) + ylab(strivar[2])
  }
  # Return the ggplot object
  return(p)
}
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
		# If discrete, coerce to character, convert to factor, replace, relevel.
		if( is.discrete(DF[, map.var]) ){
			temp.vec <- as(DF[, map.var], "character")
			temp.vec[is.na(temp.vec)] <- id.type.rp
			DF[, map.var] <- relevel(factor(temp.vec), id.type.rp)
		}
	}
	return(DF)
}
################################################################################
#' Subset points from an ordination-derived ggplot
#'
#' Easily retrieve a plot-derived \code{data.frame} with a subset of points
#' according to a threshold and method. The meaning of the threshold depends
#' upon the method. See argument description below.
#' There are many useful examples of phyloseq ordination graphics in the
#' \href{http://joey711.github.io/phyloseq/subset_ord_plot-examples}{phyloseq online tutorials}.
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
#'  \href{http://joey711.github.io/phyloseq/subset_ord_plot-examples}{phyloseq online tutorial} for this function.
#'
#'  \code{\link{plot_ordination}}
#'
#' @import ggplot2
#' @export
#' @examples 
#' ## See the online tutorials.
#' ## http://joey711.github.io/phyloseq/subset_ord_plot-examples
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
#'  \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{phyloseq online tutorials}
#'
#' @import ggplot2
#' @export
#' @examples
#' # First load and trim a dataset
#' data("GlobalPatterns")
#' GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]), GlobalPatterns)
#' # Test plots (preforms ordination in-line, then makes scree plot)
#' plot_scree(ordinate(GP, "DPCoA", "bray"))
#' plot_scree(ordinate(GP, "PCoA", "bray"))
#' # Empty return with message
#' plot_scree(ordinate(GP, "NMDS", "bray"))
#' # Constrained ordinations
#' plot_scree(ordinate(GP, "CCA", formula=~SampleType))
#' plot_scree(ordinate(GP, "RDA", formula=~SampleType)) 
#' plot_scree(ordinate(GP, "CAP", formula=~SampleType)) 
#' # Deprecated example of constrained ordination (emits a warning)
#' #plot_scree(ordinate(GP ~ SampleType, "RDA")) 
#' plot_scree(ordinate(GP, "DCA"))
#' plot_ordination(GP, ordinate(GP, "DCA"), type="scree")
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
		p = p + theme(axis.text.x=element_text(angle=90, vjust=0.5))
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
#' (instances of the phyloseq class), usually for producing graphics
#' with \code{\link{ggplot2}}. \code{psmelt} relies heavily on the 
#' \code{\link[reshape2]{melt}} and \code{\link{merge}} functions.
#' The naming conventions used in downstream phyloseq graphics functions
#' have reserved the following variable names that should not be used
#' as the names of \code{\link{sample_variables}}
#' or taxonomic \code{\link{rank_names}}.
#' These reserved names are \code{c("Sample", "Abundance", "OTU")}.
#' Also, you should not have identical names for 
#' sample variables and taxonomic ranks.
#' That is, the intersection of the output of the following two functions
#' \code{\link{sample_variables}}, \code{\link{rank_names}}
#' should be an empty vector
#' (e.g. \code{intersect(sample_variables(physeq), rank_names(physeq))}).
#' All of these potential name collisions are checked-for
#' and renamed automtically with a warning. 
#' However, if you (re)name your variables accordingly ahead of time,
#' it will reduce confusion and eliminate the warnings.
#' 
#' Note that
#' ``melted'' phyloseq data is stored much less efficiently,
#' and so RAM storage issues could arise with a smaller dataset
#' (smaller number of samples/OTUs/variables) than one might otherwise expect.
#' For common sizes of graphics-ready datasets, however,
#' this should not be a problem.
#' Because the number of OTU entries has a large effect on the RAM requirement,
#' methods to reduce the number of separate OTU entries -- 
#' for instance by agglomerating OTUs based on phylogenetic distance
#' using \code{\link{tipglom}} --
#' can help alleviate RAM usage problems.
#' This function is made user-accessible for flexibility,
#' but is also used extensively by plot functions in phyloseq.
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
#'  \code{\link[reshape2]{melt}}
#'
#'  \code{\link{merge}}
#' 
#' @import reshape2
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
  # Access covariate names from object, if present
  if(!inherits(physeq, "phyloseq")){
    rankNames = NULL
    sampleVars = NULL
  } else {
    # Still might be NULL, but attempt access
    rankNames = rank_names(physeq, FALSE)
    sampleVars = sample_variables(physeq, FALSE) 
  }
  # Define reserved names
  reservedVarnames = c("Sample", "Abundance", "OTU")  
  # type-1a conflict: between sample_data 
  # and reserved psmelt variable names
  type1aconflict = intersect(reservedVarnames, sampleVars)
  if(length(type1aconflict) > 0){
    wh1a = which(sampleVars %in% type1aconflict)
    new1a = paste0("sample_", sampleVars[wh1a])
    # First warn about the change
    warning("The sample variables: \n",
            paste(sampleVars[wh1a], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1a, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the sample variables.
    colnames(sample_data(physeq))[wh1a] <- new1a
  }
  # type-1b conflict: between tax_table
  # and reserved psmelt variable names
  type1bconflict = intersect(reservedVarnames, rankNames)
  if(length(type1bconflict) > 0){
    wh1b = which(rankNames %in% type1bconflict)
    new1b = paste0("taxa_", rankNames[wh1b])
    # First warn about the change
    warning("The rank names: \n",
            paste(rankNames[wh1b], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1b, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the conflicting taxonomic ranks
    colnames(tax_table(physeq))[wh1b] <- new1b
  }
  # type-2 conflict: internal between tax_table and sample_data
  type2conflict = intersect(sampleVars, rankNames)
  if(length(type2conflict) > 0){
    wh2 = which(sampleVars %in% type2conflict)
    new2 = paste0("sample_", sampleVars[wh2])
    # First warn about the change
    warning("The sample variables: \n",
            paste0(sampleVars[wh2], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new2, collapse=", "), "\n",
            "to avoid conflicts with taxonomic rank names.")
    # Rename the sample variables
    colnames(sample_data(physeq))[wh2] <- new2
  }
  # Enforce OTU table orientation. Redundant-looking step
  # supports "naked" otu_table as `physeq` input.
  otutab = otu_table(physeq)
  if(!taxa_are_rows(otutab)){otutab <- t(otutab)}
  # Melt the OTU table: wide form to long form table
  mdf = melt(as(otutab, "matrix"), value.name="Abundance")
  colnames(mdf)[1] <- "OTU"
  colnames(mdf)[2] <- "Sample"
  # Row and Col names are coerced to integer or factor if possible.
  # Do not want this. Coerce these to character.
  # e.g. `OTU` should always be discrete, even if OTU ID values can be coerced to integer
  mdf$OTU <- as.character(mdf$OTU)
  mdf$Sample <- as.character(mdf$Sample)
  # Merge the sample data.frame if present
  if(!is.null(sampleVars)){
    sdf = data.frame(sample_data(physeq), stringsAsFactors=FALSE)
    sdf$Sample <- sample_names(physeq)
    # merge the sample-data and the melted otu table
    mdf <- merge(mdf, sdf, by.x="Sample")
  }
  # Next merge taxonomy data, if present
  if(!is.null(rankNames)){
    TT = access(physeq, "tax_table")
    # Remove any empty columns (all NA)
    TT <- TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
    # Now add to the "psmelt" output data.frame, `mdf`
    tdf = data.frame(TT, OTU=taxa_names(physeq), stringsAsFactors=FALSE)
    mdf <- merge(mdf, tdf, by.x="OTU")
  }
  # Sort the entries by abundance
  mdf = mdf[order(mdf$Abundance, decreasing=TRUE), ]
  return(mdf)
}
################################################################################
################################################################################
#' A flexible, informative barplot phyloseq data
#'
#' There are many useful examples of phyloseq barplot graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
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
#'  \href{http://joey711.github.io/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
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
################################################################################
################################################################################
#' Returns a data table defining the line segments of a phylogenetic tree.
#'
#' This function takes a \code{\link{phylo}} or \code{\link{phyloseq-class}} object
#' and returns a list of two \code{\link{data.table}}s suitable for plotting
#' a phylogenetic tree with \code{\link[ggplot2]{ggplot}}2.
#' 
#' @param phy (Required). The \code{\link{phylo}} or \code{\link{phyloseq-class}}
#'  object (which must contain a \code{\link{phylo}}genetic tree)
#'  that you want to converted to \code{\link{data.table}}s
#'  suitable for plotting with \code{\link[ggplot2]{ggplot}}2.
#'
#' @param ladderize (Optional). Boolean or character string (either
#'  \code{FALSE}, \code{TRUE}, or \code{"left"}).
#'  Default is \code{FALSE} (no ladderization).
#'  This parameter specifies whether or not to \code{\link[ape]{ladderize}} the tree 
#'  (i.e., reorder nodes according to the depth of their enclosed
#'  subtrees) prior to plotting.
#'  This tends to make trees more aesthetically pleasing and legible in
#'  a graphical display.
#'  When \code{TRUE} or \code{"right"}, ``right'' ladderization is used.
#'  When set to \code{FALSE}, no ladderization is applied.
#'  When set to \code{"left"}, the reverse direction
#'  (``left'' ladderization) is applied.
#'  
#' @return
#'  A list of two \code{\link{data.table}}s, containing respectively 
#'  a \code{data.table} of edge segment coordinates, named \code{edgeDT},
#'  and a \code{data.table} of vertical connecting segments, named \code{vertDT}.
#'  See \code{example} below for a simple demonstration.
#' 
#' @seealso
#' An early example of this functionality was borrowed directly, with permission,
#' from the package called \code{ggphylo}, 
#' released on GitHub at:
#' \url{https://github.com/gjuggler/ggphylo}
#' by its author Gregory Jordan \email{gjuggler@@gmail.com}.
#' That original phyloseq internal function, \code{tree.layout}, has been
#' completely replaced by this smaller and much faster user-accessible 
#' function that utilizes performance enhancements from standard 
#' \code{\link{data.table}} magic as well as \code{\link{ape-package}}
#' internal C code.
#' 
#' @importFrom ape ladderize
#' @importFrom ape reorder.phylo
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @export
#' @examples
#' library("ggplot2")
#' data("esophagus")
#' phy = phy_tree(esophagus)
#' phy <- ape::root(phy, "65_2_5", resolve.root=TRUE)
#' treeSegs0 = tree_layout(phy)
#' treeSegs1 = tree_layout(esophagus)
#' edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
#' vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
#' p0 = ggplot(treeSegs0$edgeDT, edgeMap) + geom_segment() + geom_segment(vertMap, data=treeSegs0$vertDT)
#' p1 = ggplot(treeSegs1$edgeDT, edgeMap) + geom_segment() + geom_segment(vertMap, data=treeSegs1$vertDT)
#' print(p0)
#' print(p1)
#' plot_tree(esophagus, "treeonly")
#' plot_tree(esophagus, "treeonly", ladderize="left")
tree_layout = function(phy, ladderize=FALSE){
  if(inherits(phy, "phyloseq")){
    phy = phy_tree(phy)
  }
  if(!inherits(phy, "phylo")){
    stop("tree missing or invalid. Please check `phy` argument and try again.")
  }
  if(is.null(phy$edge.length)){
    # If no edge lengths, set them all to value of 1 (dendrogram).
    phy$edge.length <- rep(1L, times=nrow(phy$edge))
  }
  # Perform ladderizing, if requested
  if(ladderize != FALSE){
    if(ladderize == "left"){
      phy <- ladderize(phy, FALSE)
    } else if(ladderize==TRUE | ladderize=="right"){
      phy <- ladderize(phy, TRUE)
    } else {
      stop("You did not specify a supported option for argument `ladderize`.")
    }
  }
  # 'z' is the tree in postorder order used in calls to .C
  # Descending order of left-hand side of edge (the ancestor to the node)
  z = reorder.phylo(phy, order="postorder")
  # Initialize some characteristics of the tree.
  Nedge = nrow(phy$edge)[1]
  Nnode = phy$Nnode
  Ntip = length(phy$tip.label)
  ROOT = Ntip + 1
  TIPS = phy$edge[(phy$edge[, 2] <= Ntip), 2]
  NODES = (ROOT):(Ntip + Nnode)
  nodelabels = phy$node.label
  # Call phyloseq-internal function that in-turn calls ape's internal
  # horizontal position function, in C, using the re-ordered phylo object.
  xx = ape_node_depth_edge_length(Ntip, Nnode, z$edge, Nedge, z$edge.length)
  # Initialize `yy`, before passing to ape internal function in C.
  yy <- numeric(Ntip + Nnode)
  yy[TIPS] <- 1:Ntip
  # Define the ape_node_height wrapping function
  ape_node_height <- function(Ntip, Nnode, edge, Nedge, yy){
    .C(ape:::node_height, PACKAGE="ape",
       as.integer(Ntip), as.integer(Nnode),
       as.integer(edge[, 1]), as.integer(edge[, 2]),
       as.integer(Nedge), as.double(yy))[[6]]
  }
  # The call in ape
  #yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
  yy <- ape_node_height(Ntip, Nnode, z$edge, Nedge, yy)
  # Initialize an edge data.table 
  # Don't set key, order matters
  edgeDT = data.table(phy$edge, edge.length=phy$edge.length, OTU=NA_character_)
  # Add tip.labels if present
  if(!is.null(phy$tip.label)){
    # Initialize OTU, set node (V2) as key, assign taxa_names as OTU label
    edgeDT[, OTU:=NA_character_]
    setkey(edgeDT, V2)
    edgeDT[V2 <= Ntip, OTU:=phy$tip.label]
  }
  # Add the mapping for each edge defined in `xx` and `yy` 
  edgeDT[, xleft:=xx[V1]]
  edgeDT[, xright:=xx[V2]]
  edgeDT[, y:=yy[V2]]
  # Next define vertical segments
  vertDT = edgeDT[, list(x=xleft[1], vmin=min(y), vmax=max(y)), by=V1, mult="last"]
  if(!is.null(phy$node.label)){
    # Add non-root node labels to edgeDT
    edgeDT[V2 > ROOT, x:=xright]
    edgeDT[V2 > ROOT, label:=phy$node.label[-1]]
    # Add root label (first node label) to vertDT
    setkey(vertDT, V1)
    vertDT[J(ROOT), y:=mean(c(vmin, vmax))]
    vertDT[J(ROOT), label:=phy$node.label[1]]
  }
  return(list(edgeDT=edgeDT, vertDT=vertDT))
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
# Return TRUE if the nodes of the tree in the phyloseq object provided are unlabeled.
#' @keywords internal
nodesnotlabeled = function(physeq){
	if(is.null(phy_tree(physeq, FALSE))){
		warning("There is no phylogenetic tree in the object you have provided. Try `phy_tree(physeq)` to see.")
		return(TRUE)
	} else {
		return(is.null(phy_tree(physeq)$node.label) | length(phy_tree(physeq)$node.label)==0L)
	}
}
# A quick test function to decide how nodes should be labeled by default, if at all.
#  
#' @keywords internal
howtolabnodes = function(physeq){
	if(!nodesnotlabeled(physeq)){
    # If the nodes are labeled, use a version of this function, taking into account `ntaxa`.
		return(nodeplotdefault(manytextsize(ntaxa(physeq))))
	} else {
    # Else, use `nodeplotblank`, which returns the ggplot object as-is.
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
			p = p + geom_point(mapping=aes(x=x, y=y), data=boottop, na.rm=TRUE)
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
		p = p + geom_text(mapping=aes(x=x, y=y, label=label), data=nodelabdf,
                      size=size, hjust=hjust, na.rm=TRUE)
		return(p)
	}
}
################################################################################
#' Plot a phylogenetic tree with optional annotations
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
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
#' This function received an early development contribution from the work of 
#' Gregory Jordan via \href{https://github.com/gjuggler/ggphylo}{the ggphylo package}.
#' \code{plot_tree} has since been re-written.
#' For details see \code{\link{tree_layout}}.
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
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
#' 
#' @param shape (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point shape.
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
#'
#' @param size (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point size.
#'  A special argument \code{"abundance"} is reserved here and scales
#'  point size using abundance in each sample on a log scale.
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
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
#'  \code{FALSE}, \code{TRUE}, or \code{"left"}).
#'  Default is \code{FALSE}.
#'  This parameter specifies whether or not to \code{\link[ape]{ladderize}} the tree 
#'  (i.e., reorder nodes according to the depth of their enclosed
#'  subtrees) prior to plotting.
#'  This tends to make trees more aesthetically pleasing and legible in
#'  a graphical display.
#'  When \code{TRUE} or \code{"right"}, ``right'' ladderization is used.
#'  When set to \code{FALSE}, no ladderization is applied.
#'  When set to \code{"left"}, the reverse direction
#'  (``left'' ladderization) is applied.
#'  This argument is passed on to \code{\link{tree_layout}}.
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
#' @param treetheme (Optional).
#'  A custom \code{\link{ggplot}}2 \code{\link[ggplot2]{theme}} layer
#'  to use for the tree. Supplants any default theme layers 
#'  used within the function.
#'  A value of \code{NULL} uses a default, minimal-annotations theme. 
#'  If anything other than a them or \code{NULL}, the current global ggplot2
#'  theme will result.
#'  
#' @param justify (Optional). A character string indicating the
#'  type of justification to use on dodged points and tip labels. 
#'  A value of \code{"jagged"}, the default, results in 
#'  these tip-mapped elements being spaced as close to the tips as possible
#'  without gaps. 
#'  Currently, any other value for \code{justify} results in
#'  a left-justified arrangement of both labels and points.
#'
#' @return A \code{\link{ggplot}}2 plot.
#' 
#' @seealso
#'  \code{\link{plot.phylo}}
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
#'
#' @import scales
#' @import ggplot2
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @export
#' @examples
#' # # Using plot_tree() with the esophagus dataset.
#' # # Please note that many more interesting examples are shown
#' # # in the online tutorials"
#' # # http://joey711.github.io/phyloseq/plot_tree-examples
#' data(esophagus)
#' # plot_tree(esophagus)
#' # plot_tree(esophagus, color="Sample")
#' # plot_tree(esophagus, size="Abundance")
#' # plot_tree(esophagus, size="Abundance", color="samples")
#' plot_tree(esophagus, size="Abundance", color="Sample", base.spacing=0.03)
################################################################################
#' plot_tree(esophagus, size="abundance", color="samples", base.spacing=0.03)
plot_tree = function(physeq, method="sampledodge", nodelabf=NULL,
                       color=NULL, shape=NULL, size=NULL,
                       min.abundance=Inf, label.tips=NULL, text.size=NULL,
                       sizebase=5, base.spacing = 0.02,
                       ladderize=FALSE, plot.margin=0.2, title=NULL,
                       treetheme=NULL, justify="jagged"){
  ########################################
  # Support mis-capitalization of reserved variable names in color, shape, size
  # This helps, for instance, with backward-compatibility where "abundance"
  # was the reserved variable name for mapping OTU abundance entries
  fix_reserved_vars = function(aesvar){
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", aesvar, ignore.case=TRUE)
    return(aesvar)
  }
  if(!is.null(label.tips)){label.tips <- fix_reserved_vars(label.tips)}
  if(!is.null(color)){color <- fix_reserved_vars(color)}
  if(!is.null(shape)){shape <- fix_reserved_vars(shape)}
  if(!is.null(size) ){size  <- fix_reserved_vars(size)} 
  ########################################
  if( is.null(phy_tree(physeq, FALSE)) ){
    stop("There is no phylogenetic tree in the object you have provided.\n",
         "Try phy_tree(physeq) to see for yourself.")
  }
  if(!inherits(physeq, "phyloseq")){
    # If only a phylogenetic tree, then only tree available to overlay.
    method <- "treeonly"
  }
  # Create the tree data.table
  treeSegs <- tree_layout(phy_tree(physeq), ladderize=ladderize)
  edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
  vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
  # Initialize phylogenetic tree.
  # Naked, lines-only, unannotated tree as first layers. Edge (horiz) first, then vertical.
  p = ggplot(data=treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data=treeSegs$vertDT)
  # If no text.size given, calculate it from number of tips ("species", aka taxa)
  # This is very fast. No need to worry about whether text is printed or not.
  if(is.null(text.size)){
    text.size <- manytextsize(ntaxa(physeq))
  }
  # Add the species labels to the right.
  if(!is.null(label.tips) & method!="sampledodge"){
    # If method is sampledodge, then labels are added to the right of points, later.
    # Add labels layer to plotting object.
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if(!is.null(tax_table(object=physeq, errorIfNULL=FALSE))){
      # If there is a taxonomy available, merge it with the label data.table
      taxDT = data.table(tax_table(physeq), OTU=taxa_names(physeq), key="OTU")
      # Merge with taxonomy.
      labelDT = merge(x=labelDT, y=taxDT, by="OTU")
    }
    if(justify=="jagged"){
      # Tip label aesthetic mapping.
      # Aesthetics can be NULL, and that aesthetic gets ignored.
      labelMap <- aes_string(x="xright", y="y", label=label.tips, color=color)
    } else {
      # The left-justified version of tip-labels.
      labelMap <- aes_string(x="max(xright, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    p <- p + geom_text(labelMap, data=labelDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  }
  # Node label section.
  # 
  # If no nodelabf ("node label function") given, ask internal function to pick one.
  # Is NULL by default, meaning will dispatch to `howtolabnodes` to select function.
  # For no node labels, the "dummy" function `nodeplotblank` will return tree plot 
  # object, p, as-is, unmodified.
  if(is.null(nodelabf)){
    nodelabf = howtolabnodes(physeq)
  }
  #### set node `y` as the mean of the vertical segment
  # Use the provided/inferred node label function to add the node labels layer(s)
  # Non-root nodes first
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  # Add root label (if present)
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  # Theme specification
  if(is.null(treetheme)){
    # If NULL, then use the default tree theme.
    treetheme <- theme(axis.ticks = element_blank(),
                       axis.title.x=element_blank(), axis.text.x=element_blank(),
                       axis.title.y=element_blank(), axis.text.y=element_blank(),
                       panel.background = element_blank(),
                       panel.grid.minor = element_blank(),      
                       panel.grid.major = element_blank())   
  }
  if(inherits(treetheme, "theme")){
    # If a theme, add theme layer to plot. 
    # For all other cases, skip this, which will cause default theme to be used
    p <- p + treetheme
  }
  # Optionally add a title to the plot
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }  
  if(method!="sampledodge"){
    # If anything but a sampledodge tree, return now without further decorations.
    return(p)
  }
  ########################################
  # Sample Dodge Section
  # Special words, c("Sample", "Abundance", "OTU")
  # See psmelt()
  ########################################
  # Initialize the species/taxa/OTU data.table
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  # Merge with psmelt() result, to make all co-variables available
  dodgeDT = merge(x=dodgeDT, y=data.table(psmelt(physeq), key="OTU"), by="OTU")
  if(justify=="jagged"){
    # Remove 0 Abundance value entries now, not later, for jagged.
    dodgeDT <- dodgeDT[Abundance > 0, ]    
  }
  # Set key. Changes `dodgeDT` in place. OTU is first key, always.
  if( !is.null(color) | !is.null(shape) | !is.null(size) ){
    # If color, shape, or size is chosen, setkey by those as well
    setkeyv(dodgeDT, cols=c("OTU", color, shape, size))
  } else {
    # Else, set key by OTU and sample name. 
    setkey(dodgeDT, OTU, Sample)
  }
  # Add sample-dodge horizontal adjustment index. In-place data.table assignment
  dodgeDT[, h.adj.index := 1:length(xright), by=OTU]
  # `base.spacing` is a user-input parameter.
  # The sampledodge step size is based on this and the max `x` value
  if(justify=="jagged"){
    dodgeDT[, xdodge:=(xright + h.adj.index * base.spacing * max(xright, na.rm=TRUE))]
  } else {
    # Left-justified version, `xdodge` always starts at the max of all `xright` values.
    dodgeDT[, xdodge := max(xright, na.rm=TRUE) + h.adj.index * base.spacing * max(xright, na.rm=TRUE)]
    # zeroes removed only after all sample points have been mapped.
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  # The general tip-point map. Objects can be NULL, and that aesthetic gets ignored.
  dodgeMap <- aes_string(x="xdodge", y="y", color=color, fill=color,
                        shape=shape, size=size)
  p <- p + geom_point(dodgeMap, data=dodgeDT, na.rm=TRUE)
  # Adjust point size transform
  if( !is.null(size) ){
    p <- p + scale_size_continuous(trans=log_trans(sizebase))
  }  
  # Optionally-add abundance value label to each point.
  # User controls this by the `min.abundance` parameter.
  # A value of `Inf` implies no labels.
  if( any(dodgeDT$Abundance >= min.abundance[1]) ){
    pointlabdf = dodgeDT[Abundance>=min.abundance[1], ]
    p <- p + geom_text(mapping=aes(xdodge, y, label=Abundance),
                       data=pointlabdf, size=text.size, na.rm=TRUE)
  }
  # If indicated, add the species labels to the right of dodged points.
  if(!is.null(label.tips)){
    # `tiplabDT` has only one row per tip, the farthest horizontal
    # adjusted position (one for each taxa)
    tiplabDT = dodgeDT
    tiplabDT[, xfartiplab:=max(xdodge), by=OTU]
    tiplabDT <- tiplabDT[h.adj.index==1, .SD, by=OTU]
    if(!is.null(color)){
      if(color %in% sample_variables(physeq, errorIfNULL=FALSE)){
        color <- NULL
      }
    }
    labelMap <- NULL
    if(justify=="jagged"){
      labelMap <- aes_string(x="xfartiplab", y="y", label=label.tips, color=color)
    } else {
      labelMap <- aes_string(x="max(xfartiplab, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    # Add labels layer to plotting object.
    p <- p + geom_text(labelMap, tiplabDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  } 
  # Plot margins. 
  # Adjust the tree graphic plot margins.
  # Helps to manually ensure that graphic elements aren't clipped,
  # especially when there are long tip labels.
  min.x <- -0.01 # + dodgeDT[, min(c(xleft))]
  max.x <- dodgeDT[, max(xright, na.rm=TRUE)]
  if("xdodge" %in% names(dodgeDT)){
    max.x <- dodgeDT[, max(xright, xdodge, na.rm=TRUE)]
  }
  if(plot.margin > 0){
    max.x <- max.x * (1.0 + plot.margin)
  } 
  p <- p + scale_x_continuous(limits=c(min.x, max.x))  
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
#' \href{http://joey711.github.io/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
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
#' @param sample.order (Optional). Default \code{NULL}. 
#'  Either a single character string matching 
#'  one of the \code{\link{sample_variables}} in your data,
#'  or a character vector of \code{\link{sample_names}}
#'  in the precise order that you want them displayed in the heatmap.
#'  This overrides any ordination ordering that might be done
#'  with the \code{method}/\code{distance} arguments.
#'  
#' @param taxa.order (Optional). Default \code{NULL}. 
#'  Either a single character string matching 
#'  one of the \code{\link{rank_names}} in your data,
#'  or a character vector of \code{\link{taxa_names}}
#'  in the precise order that you want them displayed in the heatmap.
#'  This overrides any ordination ordering that might be done
#'  with the \code{method}/\code{distance} arguments.
#' 
#' @param first.sample (Optional). Default \code{NULL}.
#'  A character string matching one of the \code{\link{sample_names}}
#'  of your input data (\code{physeq}). 
#'  It will become the left-most sample in the plot.
#'  For the ordination-based ordering (recommended),
#'  the left and right edges of the axes are adjaacent in a continuous ordering. 
#'  Therefore, the choice of starting sample is meaningless and arbitrary,
#'  but it is aesthetically poor to have the left and right edge split 
#'  a natural cluster in the data.
#'  This argument allows you to specify the left edge
#'  and thereby avoid cluster-splitting, emphasize a gradient, etc.
#'  
#' @param first.taxa (Optional). Default \code{NULL}.
#'  A character string matching one of the \code{\link{taxa_names}}
#'  of your input data (\code{physeq}). 
#'  This is equivalent to \code{first.sample} (above),
#'  but for the taxa/OTU indices, usually the vertical axis.
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
#' \href{http://joey711.github.io/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
#' 
#' @importFrom vegan scores
#' @import scales 
#' 
#' @export
#' @examples
#' data("GlobalPatterns")
#' gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
#' # FYI, the base-R function uses a non-ecological ordering scheme,
#' # but does add potentially useful hclust dendrogram to the sides...
#' gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
#' # Remove the nearly-empty samples (e.g. 10 reads or less)
#' gpac = prune_samples(sample_sums(gpac) > 50, gpac)
#' # Use DESeq2 variance-stabilizing transformation of counts
#' gpacds2 = phyloseq_to_deseq2(gpac, design=~SampleType)
#' library("DESeq2")
#' gpacds2 = estimateSizeFactors(gpacds2)
#' gpacds2 = estimateDispersions(gpacds2, fitType="local")
#' gpacvst = getVarianceStabilizedData(gpacds2)
#' otu_table(gpac) <- otu_table(gpacvst, TRUE)
#' # Set values below zero, to zero.
#' otu_table(gpac)[otu_table(gpac) < 0.0] <- 0
#' # Arbitrary order if method set to NULL
#' plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family")
#' # Use ordination
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family")
#' # Use ordination for OTUs, but not sample-order
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", sample.order="SampleType")
#' # Specifying both orders omits any attempt to use ordination. The following should be the same.
#' p0 = plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", taxa.order="Phylum", sample.order="SampleType")
#' p1 = plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family", taxa.order="Phylum", sample.order="SampleType")
#' #expect_equivalent(p0, p1)
#' # Example: Order matters. Random ordering of OTU indices is difficult to interpret, even with structured sample order
#' rando = sample(taxa_names(gpac), size=ntaxa(gpac), replace=FALSE)
#' plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family", taxa.order=rando, sample.order="SampleType")
#' # # Select the edges of each axis. 
#' # First, arbitrary edge, ordering
#' plot_heatmap(gpac, method=NULL)
#' # Second, biological-ordering (instead of default ordination-ordering), but arbitrary edge
#' plot_heatmap(gpac, taxa.order="Family", sample.order="SampleType")
#' # Third, biological ordering, selected edges
#' plot_heatmap(gpac, taxa.order="Family", sample.order="SampleType", first.taxa="546313", first.sample="NP2")
#' # Fourth, add meaningful labels
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", taxa.order="Family", sample.order="SampleType", first.taxa="546313", first.sample="NP2")
plot_heatmap <- function(physeq, method="NMDS", distance="bray", 
	sample.label=NULL, taxa.label=NULL, 
	low="#000033", high="#66CCFF", na.value="black", trans=log_trans(4), 
	max.label=250, title=NULL, sample.order=NULL, taxa.order=NULL,
  first.sample=NULL, first.taxa=NULL, ...){

  # User-override ordering
  if( !is.null(taxa.order) & length(taxa.order)==1 ){
    # Assume `taxa.order` is a tax_table variable. Use it for ordering.
    rankcol = which(rank_names(physeq) %in% taxa.order)
    taxmat = as(tax_table(physeq)[, 1:rankcol], "matrix")
    taxa.order = apply(taxmat, 1, paste, sep="", collapse="")
    names(taxa.order) <- taxa_names(physeq)
    taxa.order = names(sort(taxa.order, na.last=TRUE))
  }
  if( !is.null(sample.order) & length(sample.order)==1 ){
    # Assume `sample.order` is a sample variable. Use it for ordering.
    sample.order = as.character(get_variable(physeq, sample.order))
    names(sample.order) <- sample_names(physeq)
    sample.order = names(sort(sample.order, na.last=TRUE))
  }
  
  if( !is.null(method) & (is.null(taxa.order) | is.null(sample.order)) ){
    # Only attempt NeatMap if method is non-NULL & at least one of
    # taxa.order and sample.order is not-yet defined.
    # If both axes orders pre-defined by user, no need to perform ordination...
    
	  # Copy the approach from NeatMap by doing ordination on samples, but use 
	  # phyloseq-wrapped distance/ordination procedures.
	  # Reorder by the angle in radial coordinates on the 2-axis plane.
    
    # Capture the NMDS iterations cat() output with capture.output
		junk = capture.output( ps.ord <- ordinate(physeq, method, distance, ...), file=NULL)
    if( is.null(sample.order) ){
      # Only define new ord-based sample order if user did not define one already
      reduction.result = scores(ps.ord, choices=c(1, 2), display="sites")
      sample.order = sample_names(physeq)[order(RadialTheta(reduction.result))]      
    }

		test <- try(scores(ps.ord, choices=c(1, 2), display="species"), TRUE)
		if( class(test) != "try-error" & !is.null(test) & is.null(taxa.order) ){			
		  # re-order species/taxa/OTUs, if possible,
		  # and only if user did not define an order already
			OTUreduct = scores(ps.ord, choices=c(1, 2), display="species")
			taxa.order  = taxa_names(physeq)[order(RadialTheta(OTUreduct))]
		}
	}
  
  # Now that index orders are determined, check/assign edges of axes, if specified
  if( !is.null(first.sample) ){
    sample.order = restart(sample.order, first.sample)
  }
  if( !is.null(first.taxa) ){
    taxa.order = restart(taxa.order, first.taxa)
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
	if( !is.null(taxa.order) ){
		# If OTU-order is available, coerce to factor with special level-order
		adf$OTU = factor(adf$OTU, levels=taxa.order)
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
	  # Remove the labels from any rendering.
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
		if( !is.null(taxa.order) ){		
			# Re-order according to taxa.order
			labvec <- labvec[taxa.order]
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
# Chunk re-order a vector so that specified newstart is first.
# Different than relevel.
#' @keywords internal
restart = function(x, newstart){
  # x = sample_names(gpac)
  # newstart = "NP2"
  pivot = which(x %in% newstart)
  return(c(x[pivot:length(x)], x[1:(pivot-1)]))
}
################################################################################
#' Create a ggplot summary of gap statistic results
#'
#' @param clusgap (Required). 
#' An object of S3 class \code{"clusGap"}, basically a list with components.
#' See the \code{\link[cluster]{clusGap}} documentation for more details.
#' In most cases this will be the output of \code{\link{gapstat_ord}},
#' or \code{\link[cluster]{clusGap}} if you called it directly.
#' 
#' @param title (Optional). Character string.
#'  The main title for the graphic.
#'  Default is \code{"Gap Statistic results"}.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} plot object. 
#' The rendered graphic should be a plot of the gap statistic score 
#' versus values for \code{k}, the number of clusters.
#' 
#' @seealso
#' \code{\link{gapstat_ord}}
#' 
#' \code{\link[cluster]{clusGap}}
#' 
#' \code{\link[ggplot2]{ggplot}}
#' 
#' @import ggplot2
#' @export
#' @examples
#' # Load and process data
#' data("soilrep")
#' soilr = rarefy_even_depth(soilrep, rngseed=888)
#' print(soilr)
#' sample_variables(soilr)
#' # Ordination
#' sord  = ordinate(soilr, "DCA")
#' # Gap Statistic
#' gs = gapstat_ord(sord, axes=1:4, verbose=FALSE)
#' # Evaluate results with plots, etc.
#' plot_scree(sord)
#' plot_ordination(soilr, sord,  color="Treatment")
#' plot_clusgap(gs)
#' print(gs, method="Tibs2001SEmax")
#' # Non-ordination example, use cluster::clusGap function directly
#' library("cluster")
#' pam1 = function(x, k){list(cluster = pam(x, k, cluster.only=TRUE))}
#' gs.pam.RU = clusGap(ruspini, FUN = pam1, K.max = 8, B = 60)
#' gs.pam.RU
#' plot(gs.pam.RU, main = "Gap statistic for the 'ruspini' data")
#' mtext("k = 4 is best .. and  k = 5  pretty close")
#' plot_clusgap(gs.pam.RU)
plot_clusgap = function(clusgap, title="Gap Statistic results"){
	gstab = data.frame(clusgap$Tab, k = 1:nrow(clusgap$Tab))
	p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
	p = p + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim))
	p = p + ggtitle(title)
	return(p)
}
################################################################################