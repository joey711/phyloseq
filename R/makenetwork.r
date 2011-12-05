###############################################################################
#' Create a taxa graph adjacency matrix from an \code{otuTable}.
#' 
#' Function that inputs an abundance table as a matrix 
#' and makes a graph adjacency matrix ready to plot 
#' with the igraph package.
#'
#' @usage makenetwork(abund, plotgraph=TRUE, 
#'			community=TRUE, threshold=0, incommon=0.4, method="jaccard")
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
#'
#' @export
#' @rdname makenetwork-methods
#' @docType methods
#'
#' @examples #
#' ## data(ex1)
#' ## makenetwork(otuTable(ex1), TRUE)
setGeneric("makenetwork", function(abund, plotgraph=TRUE, 
			community=TRUE, threshold=0, incommon=0.4, method="jaccard"){
	standardGeneric("makenetwork")				
})
###############################################################################
# Creates a taxa network from an otuTable.
#' @import igraph vegan
#' @aliases makenetwork,otuTable-method
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
	jaccpa <- vegdist(presenceAbsence, method)
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
# Creates a taxa network from an otuTable.
#' @aliases makenetwork,phyloseq-method
#' @rdname makenetwork-methods
setMethod("makenetwork", signature("phyloseq"),
	function(abund,plotgraph=TRUE,community=TRUE,threshold=0,incommon=0.4,method="jaccard"){
		makenetwork(t(otuTable(abund)), plotgraph, community, threshold, incommon, method)
})
###############################################################################
################################################################################
#' Convert edgelist hash-table to clique list
#'
#' Agglomerating function to convert an edgelist -- which is a 2-column table
#' of vertices where each row represents an edge -- to a list of cliques,
#' in which each clique is represented by a character vector of the vertex labels
#' of the vertices that are members of the clique. This algorithm is perfectly
#' greedy, such that the only requirement for inclusion in a clique is an edge
#' to any of the other members of that clique.
#'
#' @usage edgelist2clique(EdgeList)
#'
#' @param EdgeList a 2-column table of vertices where each row represents an edge. 
#'
#' @return A list, where each element is a character vector of tips that should
#' are in the same clique.
#'
#' @export
#' @examples #
#' # edgelist2clique(get.edgelist(ig))
edgelist2clique = function(EdgeList){
	# initialize list of globs
	glob	= vector(mode="list")
	
	for (i in 1:nrow(EdgeList) ){
		# initialize for each loop a 'skip' variable, to avoid later tests if answer already found.
		skip	= FALSE
		# check entries in list, glob, for membership in already-forming glob
		thislink	= unlist(EdgeList[i,1:2])
		# identify which, if any, globs have contig of interest already in them
		# OLD WAY 1
		# glob1	= which(sapply(sapply(glob,function(globi,ctig){which(ctig==globi)},ctig=thislink[1]),length)>0)
		# glob2	= which(sapply(sapply(glob,function(globi,ctig){which(ctig==globi)},ctig=thislink[2]),length)>0)
		# better way
		glob1 = which(vapply(glob, function(glb,vertex){vertex %in% glb}, TRUE, thislink[1]))
		glob2 = which(vapply(glob, function(glb,vertex){vertex %in% glb}, TRUE, thislink[2]))
		# grep is actually slower
		# glob1 = grep(thislink[1], glob, fixed=TRUE)
		# glob2 = grep(thislink[2], glob, fixed=TRUE)
	
		## Now series of if-tests to decide where and how to glob.
		# if both contigs are already in same globs, skip
		if( !skip & length(glob1)>0 & length(glob2)>0 ){
			if( glob1==glob2 ){
				# skip remaining tests, only.
				skip = TRUE	
			}
		}
		# if both contigs are already in different globs, join both globs
		if( !skip & length(glob1)>0 & length(glob2)>0 ){
			if( glob1!=glob2 ){
				# join globs at end	
				glob	= c(glob,list(c(glob[[glob1]],glob[[glob2]])))
				# remove old globs
				glob	= glob[-c(glob1,glob2)]
				# skip remaining tests.
				skip = TRUE	
			}
		}
		# if only contig 1 is in a glob, add Contig 2 to that glob
		if( !skip & length(glob1)>0 & length(glob2)==0){
			# add Contig 2 to glob1	
			glob[[glob1]]	= c(glob[[glob1]],thislink[2])
			# skip remaining tests.
			skip = TRUE	
		}
		# if only Contig 2 is in a glob, add contig 1 to that glob
		if( !skip & length(glob2)>0 & length(glob1)==0 ){
			# add Contig 1 to glob2	
			glob[[glob2]]	= c(glob[[glob2]],thislink[1])
			# skip remaining tests.
			skip = TRUE		
		}
		# all else, form new glob
		if( !skip & length(glob1)==0 & length(glob2)==0 ){
			# add Contig 1 and Contig 2 as vector in new glob.
			glob	= c(glob,list(thislink))
		}
	}
	return(glob)
}
################################################################################