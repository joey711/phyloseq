################################################################################
#' Make microbiome network (igraph0)
#'
#' A specialized function for creating a network representation of microbiomes,
#' sample-wise or taxa-wise,
#' based on a user-defined ecological distance and (potentially arbitrary) threshold.
#' The graph is ultimately represented using the 
#' \code{igraph0}-package.
#'
#' @usage make_network(physeq, type="samples", distance="jaccard", max.dist = 0.4, 
#'     keep.isolates=FALSE, ...)
#'
#' @param physeq (Required). Default \code{NULL}. 
#'  A \code{\link{phyloseq-class}} object,
#'  or \code{\link{otuTable-class}} object,
#'  on which \code{g} is based. \code{phyloseq-class} recommended.
#'
#' @param type (Optional). Default \code{"samples"}.
#'  Whether the network should be samples or taxa/species.
#'  Supported arguments are \code{"samples"}, \code{"species"},
#'  where \code{"species"} indicates using the taxa indices,
#'  whether they actually represent species or some other taxonomic rank.
#'
#'  NOTE: not all distance methods are supported if \code{"species"}
#'  selected for type. For example, the UniFrac distance and DPCoA
#'  cannot be calculated for taxa-wise distances, because they use
#'  a taxa-wise tree as part of their calculation between samples, and
#'  there is no transpose-equivalent for this tree.
#'
#' @param distance (Optional). Default \code{"jaccard"}.
#'  Any supported argument to the \code{method} parameter of the 
#'  \code{\link{distance}} function is supported here.
#'  Some distance methods, like \code{"unifrac"}, may take 
#'  a non-trivial amount of time to calculate, in which case
#'  you probably want to calculate the distance matrix separately,
#'  save, and then provide it as the argument to \code{distance} instead.
#'  See below for alternatives).
#'
#'  Alternatively, if you have already calculated the sample-wise distance
#'  object, the resulting \code{dist}-class object
#'  can be provided as \code{distance} instead (see examples).
#' 
#'  A third alternative is to provide a function that takes 
#'  a sample-by-species matrix (typical vegan orientation)
#'  and returns a sample-wise distance
#'  matrix. 
#' 
#' @param max.dist (Optional). Default \code{0.4}. 
#'  The maximum ecological distance (as defined by \code{distance})
#'  allowed between two samples to still consider them ``connected'' 
#'  by an edge in the graphical model.
#' 
#' @param keep.isolates (Optional). Default \code{FALSE}. Logical.
#'  Whether to keep isolates (un-connected samples, not microbial isolates)
#'  in the graphical model that is returned. Default results in isolates
#'  being removed from the object.
#'
#' @param ... (Optional). Additional parameters passed on to \code{\link{distance}}.
#'
#' @return A \code{igraph0}-class object. 
#' 
#' @seealso 
#'  \code{\link{plot_network}}
#'
#' @importFrom igraph0 graph.adjacency
#' @importFrom igraph0 V
#' @importFrom igraph0 delete.vertices
#' @importFrom igraph0 degree
#'
#' @export
#'
#' @examples
#' # # Example plots with Enterotype Dataset
#' data(enterotype)
#' ig <- make_network(enterotype, max.dist=0.3)
#' plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # 
#' # ig <- make_network(enterotype, max.dist=0.2)
#' # plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # 
#' # # Three methods of choosing/providing distance/distance-method
#' # Provide method name available to distance
#' ig <- make_network(enterotype, max.dist=0.3, distance="jaccard")
#' # Provide distance object, already computed
#' jaccdist <- distance(enterotype, "jaccard")
#' ih <- make_network(enterotype, max.dist=0.3, distance=jaccdist)
#' # Provide "custom" function.
#' ii <- make_network(enterotype, max.dist=0.3, distance=function(x){vegdist(x, "jaccard")})
#' # The have equal results:		
#' all.equal(ig, ih)
#' all.equal(ig, ii)
#' # 
#' # Try out making a trivial "network" of the 3-sample esophagus data,
#' # with weighted-UniFrac as distance
#' data(esophagus)
#' ij <- make_network(esophagus, "samples", "unifrac", weighted=TRUE)
make_network <- function(physeq, type="samples", distance="jaccard", max.dist = 0.4, 
	keep.isolates=FALSE, ...){

	if( type == "species"){
	    # Calculate or asign species-wise distance matrix
	    if( class(distance) == "dist" ){ # If argument is already a distance matrix.
	    	# If distance a distance object, use it rather than re-calculate
	    	obj.dist <- distance
	    	if( attributes(obj.dist)$Size != nspecies(physeq) ){
	    		stop("nspecies(physeq) does not match size of dist object in distance")
	    	}
	    	if( !setequal(attributes(obj.dist)$Labels, species.names(physeq)) ){
	    		stop("species.names does not exactly match dist-indices")
	    	}
	
	    # If character string, pass on to distance(), assume supported
	    } else if( class(distance) == "character" ){ 
		    obj.dist <- distance(physeq, method=distance, type=type, ...)
	
		# Else, assume a custom function and attempt to calculate.
	    } else { 	
	    	# Enforce orientation for species-wise distances
		    if( !speciesAreRows(physeq) ){ physeq <- t(physeq) }
		    
		    # Calculate distances
		    obj.dist <- distance(as(otuTable(physeq), "matrix"))	    
	    }
	    # coerce distance-matrix back into vanilla matrix, Taxa Distance Matrix, TaDiMa
	    TaDiMa  <- as.matrix(obj.dist)
	    
	    # Add Inf to the diagonal to avoid self-connecting edges (inefficient)
	    TaDiMa <- TaDiMa + diag(Inf, nspecies(physeq), nspecies(physeq))
	    
	    # Convert distance matrix to coincidence matrix, CoMa, using max.dist
		CoMa <- TaDiMa < max.dist   
			
	} else if( type == "samples" ){
		
	    # Calculate or asign sample-wise distance matrix
	    if( class(distance) == "dist" ){ # If argument is already a distance matrix.
	    	# If distance a distance object, use it rather than re-calculate
	    	obj.dist <- distance
	    	if( attributes(obj.dist)$Size != nsamples(physeq) ){
	    		stop("nsamples(physeq) does not match size of dist object in distance")
	    	}
	    	if( !setequal(attributes(obj.dist)$Labels, sample.names(physeq)) ){
	    		stop("sample.names does not exactly match dist-indices")
	    	}
	    	
	    # If character string, pass on to distance(), assume supported    	
	    } else if( class(distance) == "character" ){
		    obj.dist <- distance(physeq, method=distance, type=type, ...)
	
		# Else, assume a custom function and attempt to calculate.	    
	    } else { 
	    	# Enforce orientation for sample-wise distances
		    if(speciesAreRows(physeq)){ physeq <- t(physeq) }
		    
		    # Calculate distances
		    obj.dist <- distance(as(otuTable(physeq), "matrix"))	    
	    }
	    # coerce distance-matrix back into vanilla matrix, Sample Distance Matrix, SaDiMa
	    SaDiMa  <- as.matrix(obj.dist)
	    
	    # Add Inf to the diagonal to avoid self-connecting edges (inefficient)
	    SaDiMa <- SaDiMa + diag(Inf, nsamples(physeq), nsamples(physeq))
	    
	    # Convert distance matrix to coincidence matrix, CoMa, using max.dist
		CoMa <- SaDiMa < max.dist  
		  
	} else {
		stop("type argument must be one of \n (1) samples \n or \n (2) species")
	}
    
    # Calculate the igraph0-formatted network
    ig <- graph.adjacency(CoMa, mode="lower")
    
    # If keeping isolates, done. Else, remove them, then return igraph0.
    if( keep.isolates ){
    	return(ig)
    } else {
	    isolates   <- V(ig)[degree(ig) == 0]
	    ig.no.isol <- delete.vertices(ig, V(ig)[degree(ig) == 0])
	    return(ig.no.isol)
    }
}
################################################################################