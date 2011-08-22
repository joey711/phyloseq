###############################################################################
#' Create a taxa graph adjacency matrix from an \code{otuTable}.
#' 
#' Function that inputs an abundance table as a matrix 
#' and makes a graph adjacency matrix ready to plot 
#' with the igraph package.
#'
#' @param abund The species abundance table (otuTable) or a more complex
#'  phyloseq-package object that contains an otuTable.
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
#' @export
#' @rdname makenetwork-methods
#' @examples #
#' ## data(ex1)
#' ## makenetwork(otuTable(ex1), TRUE)
setGeneric("makenetwork", function(abund, plotgraph=TRUE, 
			community=TRUE, threshold=0, incommon=0.4, method="jaccard")
	standardGeneric("makenetwork")
)
###############################################################################
#' Creates a taxa network from an otuTable.
#'
#' @import igraph 
#' @author Susan Holmes
#' @exportMethod makenetwork
#' @examples #
#'
#' @aliases makenetwork,otuTable-method
#' @docType methods
#' @rdname makenetwork-methods
setMethod("makenetwork", "otuTable", function(abund, plotgraph=TRUE, 
	community=TRUE, threshold=0, incommon=0.4, method="jaccard"){
	# abundance is the abundance table with no rows that are all zero
	# plotgraph is a toggle for whether a plotted network is requested
	# if commun=TRUE the function will also output the community groups	
	# threshold is the number that is fixed for when prsence occurs
	# for instance, threshold>1 means that species with 2 or more reads
	# are considered present	
	#require(vegan); require(igraph)
	### Keep the original row numbers as labels
	### this is for the later analysis of groups
	if (is.na(dimnames(abund)[[1]][1])) dimnames(abund)=list(1:nrow(abund),1:ncol(abund))
	### Only take the rows where there are at least one value over threshold
	abundance <- abund[rowSums(abund)>threshold,]
	n         <- nrow(abundance)

	# Convert to 1,0 binary matrix for input to vegdist. -0 converts to numeric
	presenceAbsence <- (abundance > threshold) - 0

	##Compute the Jaccard distance between the rows, this will only make points
	##closer if they are actually present together	
	##You could use any of the other distances in vegan or elsewhere
	jaccpa <- vegan::vegdist(presenceAbsence, method)
	###Distances in R are vectors by default, we make them into matrices	
	jaacm <- as.matrix(jaccpa)
	coinc <- matrix(0,n,n)
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
})
###############################################################################
# Extend for PhyloSeq data class
###############################################################################
#' Creates a taxa network from an otuTable.
#'
#' @author Susan Holmes
#' @exportMethod makenetwork
#' @examples #
#'
#' @aliases makenetwork,phyloseq-method
#' @docType methods
#' @rdname makenetwork-methods
setMethod("makenetwork", signature("phyloseq"),
	function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="jaccard"){
		callNextMethod(t(otuTable(abund)),plotgraph,community,threshold,incommon,method)
})
###############################################################################
#' Creates a taxa network from an otuTable.
#'
#' @author Susan Holmes
#' @exportMethod makenetwork
#' @examples #
#'
#' @aliases makenetwork,otuTree-method
#' @docType methods
#' @rdname makenetwork-methods
setMethod("makenetwork", signature("otuTree"),
	function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="jaccard"){
		callNextMethod(t(otuTable(abund)),plotgraph,community,threshold,incommon,method)
})
###############################################################################