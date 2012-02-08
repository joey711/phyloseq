\docType{methods}
\name{taxglom}
\alias{taxglom}
\alias{taxglom,otuTable,character-method}
\alias{taxglom,otuTable,taxonomyTable-method}
\alias{taxglom,phyloseq,ANY-method}
\title{Agglomerate taxa of the same type.}
\usage{
  taxglom(physeq, tax=NULL, taxlevel="Phylum", NArm=TRUE,
  bad_empty=c(NA, "", " ", "\t"))
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}} or
  \code{\link{otuTable}}.}

  \item{tax}{(Optional). Either a
  \code{link{taxonomyTable-class}}, or alternatively, a
  character vector specifying the desired taxonomic group
  of each taxa in \code{physeq}. If \code{tax} is a
  character vector, it must have length equal to the
  (original) number of taxa in \code{physeq}
  (\code{nspecies(physeq)}), and each element must be named
  according to the taxa ID (that is, the result of
  \code{species.names(physeq)}). If \code{tax} is a
  character vector, than the \code{taxlevel} argument is
  ignored. If \code{physeq} already contains a
  \code{taxonomyTable} component in its \code{taxTab} slot,
  then the \code{tax} argument is ignored.}

  \item{taxlevel}{A single-element character specifying the
  taxonomic level (column name) in \code{tax}, the
  \code{taxonomyTable}, that you want to agglomerate over.
  The default value is \code{"Phylum"}. Note that this
  default may agglomerate too broadly for a given
  experiment, and the user is strongly encouraged to try
  different taxonomic levels.}

  \item{NArm}{(Optional). Logical, length equal to one.
  Default is \code{TRUE}. CAUTION. The decision to prune
  (or not) taxa for which you lack categorical data could
  have a large effect on downstream analysis. You may want
  to re-compute your analysis under both conditions, or at
  least think carefully about what the effect might be and
  the reasons explaining the absence of information for
  certain taxa. In the case of taxonomy, it is often a
  result of imprecision in taxonomic designation based on
  short phylogenetic sequences and a patchy system of
  nomenclature. If this seems to be an issue for your
  analysis, think about also trying the
  nomenclature-agnostic \code{\link{tipglom}} method if you
  have a phylogenetic tree available.}

  \item{bad_empty}{(Optional). Character vector. Default:
  \code{c(NA, "", " ", "\t")}. Defines the bad/empty values
  that should be ignored and/or considered unknown. They
  will be removed from the internal agglomeration vector
  derived from the argument to \code{tax}, and therefore
  agglomeration will not combine taxa according to the
  presence of these values in \code{tax}. Furthermore, the
  corresponding taxa can be optionally pruned from the
  output if \code{NArm} is set to \code{TRUE}.}
}
\value{
  A taxonomically-agglomerated, optionally-pruned, object
  with class matching the class of \code{physeq}.
}
\description{
  This method merges species if, at a certain taxaonomic
  rank, their taxonomy is the same. Its approach is
  analogous to \code{tipglom}, but uses categorical data
  instead of a tree. In principal, other categorical data
  known for all taxa could also be used in place of
  taxonomy.
}
\examples{
# data(ex1)
# ## print the available taxonomic ranks
# colnames(taxTab(ex1))
# ## agglomerate at the Family taxonomic rank
# (x1 <- taxglom(ex1, taxlevel="Family") )
# ## How many taxa before/after agglomeration?
# nspecies(ex1); nspecies(x1)
# ## Look at enterotype dataset...
# data(enterotype)
# ## print the available taxonomic ranks. Shows only 1 rank available, not useful for taxglom
# colnames(taxTab(enterotype))
}
\seealso{
  \code{\link{tipglom}}, \code{\link{prune_species}},
  \code{\link{merge_species}}
}

