################################################################################
#' Wrapper for parallelization
#'
#' This function is a convenience wrapper around \code{Rmpi::mpi.spawn.Rslaves()}
#' and especially \code{Rmpi::mpi.apply}. In essence, it should act as a parallel
#' apply-loop on the iteration vector, \code{i}.
#'
#' @usage parapply(i, fun1, num_cores, ...)
#'
#' @param i (Required). The vector/list to loop through.
#' @param fun1 (Required). The function to apply.
#' @param num_cores (Required). Integer. The number of cores to use.
#' @param ... (Optional). Additional arguments, passed to fun1.
#' 
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' @keywords internal
#' 
parapply = function(i, fun1, num_cores, ...){
	require(Rmpi)
	# define number of slaves for chunking, ns:
	ns 	= num_cores - 1
	# define outerloop, OL. Has length `ns' and contains the list of indices
	# for each paral job
	e	= 1:ns*ceiling((length(i)/ns))
	s	= e - e[1] + 1
	# trim indexes so that outerloop matches index length
	e[length(e)] <- length(i) 
	# create outer loop: the set of indexes to send to fun1 via mpi.apply
	OL	= lapply(1:length(s),function(j,s,e){s[j]:e[j]},s,e)
	# now spawn your slaves.
	Rmpi::mpi.spawn.Rslaves(nslaves=ns)
	# Finally, run your chunked loop.
	para_list <- Rmpi::mpi.apply(OL, fun1, ...)
	Rmpi::mpi.close.Rslaves()
	return(unlist(para_list))
}
################################################################################