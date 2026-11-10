# The two trees differ only by a zero-length resolution of a polytomy.
# UniFrac should therefore give identical distances for both trees.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else "."
pkgload::load_all(root, quiet=TRUE)
counts <- otu_table(matrix(c(1, 0, 0, 0, 0, 0, 1, 0), nrow=4,
                           dimnames=list(c("a", "b", "c", "d"), c("u", "v"))), TRUE)
tree <- ape::read.tree(text="((a:1,b:1,c:1):1,d:2);")
resolved <- ape::multi2di(tree, random=FALSE)
raw <- phyloseq(counts, tree)
binary <- phyloseq(counts, resolved)
stopifnot(ape::is.rooted(tree), !ape::is.binary(tree),
          sum(tree$edge.length) == sum(resolved$edge.length))
collect <- function(ps, weighted, normalized=TRUE) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(as.numeric(UniFrac(ps, weighted=weighted,
                                                          normalized=normalized)),
                    warning=function(w) {
                      warnings <<- c(warnings, conditionMessage(w))
                      invokeRestart("muffleWarning")
                    }), error=function(e) conditionMessage(e))
  list(value=value, warnings=warnings)
}
result <- list(
  unweighted=list(expected=2/3, raw=collect(raw, FALSE), resolved=collect(binary, FALSE)),
  weighted_normalized=list(expected=0.5, raw=collect(raw, TRUE), resolved=collect(binary, TRUE)),
  weighted_unnormalized=list(expected=2, raw=collect(raw, TRUE, FALSE),
                            resolved=collect(binary, TRUE, FALSE)))
cat(jsonlite::toJSON(result, auto_unbox=TRUE, pretty=TRUE), "\n")
