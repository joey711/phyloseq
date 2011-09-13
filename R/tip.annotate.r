################################################################################
# tipsymbols and tiptext. Need to make both tipsymbols and tiptext documentation
# point to this one.
################################################################################
#' Annotate tips on a tree.
#'
#' There were some unexpected behavior from the tiplabels function in ape. These
#' functions are intended to act as simplified versions that act as a convenience
#' wrapper for \code{points()} or \code{text()} functions, respectively, 
#' but where the tip
#' coordinates are specified by giving the tip ID (integer) as input.
#'
#' @param tip An integer specifying the tip ID in a tree that for which the 
#'  base plot has already been generated and is still available to \code{R}.
#' 
#' @param adj A 2 element numeric vector specifying a position adjustment.
#' 
#' @param ... Additional plotting parameters that are passed to 
#'  \code{\link{points}} or \code{\link{text}} in the R base graphics.
#'
#' @return No objects returned. Symbol or text is plotted on the available
#'  graphic device.
#'
#' @export
#' @rdname tip-annotate
#' 
#' @seealso tiplabels points text
#' @examples #
#' ## data(ex1)
#' ## # for reproducibility
#' ## set.seed(711)
#' ## ex2 <- prune_species(sample(species.names(ex1), 50), ex1)
#' ## plot(as(tre(ex2), "phylo"))
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
#' Annotate the tips of a phylogenetic tree plot with text.
#'
#' Do not forget to include a \code{labels} argument.
#'
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
#' Function to plot otuSamTree with annotated tips.
#'
#' @param object A phyloseq object that contains a tree, sampleMap, and otuTable.
#'  That is, you will most-likely get an error unless this is a 
#'  otuSamTree or otuSamTaxTree object.
#'
#' @param color_factor A character string specifying the column
#'  of the sampleMap that will be used for setting the color of symbols.
#'
#' @param shape_factor A character string specifying the column
#'  of the sampleMap that will be used for setting the shape of symbols.
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
#'  be a schale, not an aesthetic map. Therefore, it will in most-cases 
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
#' @param printTheseTaxa a character vector of the taxa names in \code{object}
#'  that should be labeled on the tree plot adjacent to the right. Default is
#'  NULL. Not yet implemented.
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
plot_tree_phyloseq <- function(object, color_factor=NULL, shape_factor=NULL, 
	base_size=1, size_scaling_factor = 0.2, opacity=2/3,
	custom_color_scale=NULL, custom_shape_scale=NULL, 
	type_abundance_value=FALSE, printTheseTaxa=NULL, ...){

	# Initialize categories vector, giving the sampleMap index of aesthetics.
	categories <- character(2); names(categories) <- c("color", "shape")
	
	# Only look for default factors if color_factor AND shape_factor absent.
	if( is.null(color_factor) & is.null(shape_factor) ){
		# Want to use up to the first two factors in sampleMap as default.
		smdf_factor_types <- lapply(data.frame(sampleMap(object)), class) == "factor"
		available_factors <- colnames(sampleMap(object))[	smdf_factor_types ]
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
		color_factor <- data.frame(sampleMap(object))[, categories["color"]]
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
		color_vec <- rep("black", nsamples(object))
	}
	names(color_vec) <- sample.names(object)

	# shape	 custom_shape_scale
	if( categories["shape"] != "" ){
		shape_factor <- data.frame(sampleMap(object))[, categories["shape"]]
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
		shape_vec <- rep(22, nsamples(object))
	}
	names(shape_vec) <- sample.names(object)
	
	# This is based on ape-package plotting. Use phylo-class tree.
	tree <- as(tre(object), "phylo")
	
	# Now plot the initial, symbol-less tree. Must be first to get the proper
	# x, y limits to calculate the scales of the annotation objects.
	plot(tree, type="phylogram", show.tip.label=FALSE,
		xpd=NA, no.margin=TRUE, edge.width=1.5)

	# Store information about the tree-plot
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	#str(lastPP)
	xlims <- lastPP$x.lim
	ylims <- lastPP$y.lim	

	# Add title of tree
	treeTitle = "Tree Example"
	text(x=(0.35*xlims[2]), y=(1.02*ylims[2]), labels=treeTitle, pos=4, cex=1.5)

	# Add scale bar
	add.scale.bar(x=(0.5*xlims[2]), y=0, length=0.1)
	adj.j.start <- 0.5 + (0.01 * xlims[2])
	adj.step    <- 0.0155 * xlims[2]
	########################################
	# Now loop over sample variables
	########################################
	for( i in species.names(object) ){ # i loops over species
		#i = species.names(object)[speciesSums(object) == max(speciesSums(object))]
		# sitesi should hold the sites that are relevant to species "i"
		sitesi <- getSamples( object, i)
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
# library("phyloseq")
# data(ex1)
# ex2 <- ex1
# species.names(ex2) <- sample(species.names(ex1), 50)
# plot_tree_phyloseq(ex2)
# plot_tree_phyloseq(ex2, shape_factor="Diet")
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet")
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet", 
	# size_scaling_factor=0.6, type_abundance_value=TRUE)
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet", 
	# size_scaling_factor=0.6, custom_color_scale=c("blue", "magenta"))
