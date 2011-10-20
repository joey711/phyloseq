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
#' @keywords agglomerate clique
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