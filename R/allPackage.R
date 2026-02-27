###############################################
#'  Handling and analysis of high-throughput phylogenetic sequence data.
#'
#' There are already several ecology and phylogenetic packages available in R,
#' including the adephylo, vegan, ade4, picante, ape, phangorn, phylobase, and OTUbase packages.
#' These can already take advantage of many of the powerful statistical and graphics tools
#' available in R. However, prior to \emph{phyloseq} a user must devise their own methods
#' for parsing the output of their favorite OTU clustering application, and, as a consequence,
#' there is also no standard within Bioconductor (or R generally) for storing or sharing the
#' suite of related data objects that describe a phylogenetic sequencing project.
#' The phyloseq package seeks to address these issues by providing a related set of S4 classes
#' that internally manage the handling tasks associated with organizing, linking, storing,
#' and analyzing phylogenetic sequencing data. \emph{phyloseq} additionally provides some
#' convenience wrappers for input from common clustering applications, common analysis pipelines,
#' and native implementation of methods that are not available in other R packages.
#'
#' @import methods
#' @importFrom graphics axis
#' @importFrom stats aggregate as.dist as.formula as.hclust complete.cases cutree relevel
#' @importFrom utils capture.output combn download.file head read.table tail untar unzip write.table
#' @author Paul J. McMurdie II \email{mcmurdie@@stanford.edu}
#' @references \url{www.stanford.edu/~mcmurdie}
#' @keywords package
"_PACKAGE"

# Suppress R CMD check NOTEs for variables used in data.table NSE (:=, .SD, J),
# ggplot2 aes() column-name arguments, and other non-standard evaluation contexts.
utils::globalVariables(c(
  # data.table NSE operators and helpers
  ":=", ".SD", "J",
  # data.table / import column name variables
  "#OTU ID", "Consensus Lineage",
  "count", "queryString", "queryID", "Classification", "OTULabel",
  "read",
  # tree_layout / plot_tree column names
  "OTU", "V1", "V2", "xleft", "xright", "x", "y", "label",
  "vmin", "vmax", "h.adj.index", "xdodge", "xfartiplab",
  # plot_heatmap / plot_tree aesthetic variables
  "Sample", "Abundance",
  # plot_net / plot_network aesthetic variables
  "xend", "yend",
  # plot_richness aesthetic variables
  "value", "se",
  # plot_scree aesthetic variable
  "eigenvalue",
  # plot_clusgap aesthetic variables
  "k", "gap", "SE.sim",
  # nodeplot anonymous function variables
  # (x, y, label already listed above)
  # JSD loop variable
  "i",
  # merge_phyloseq_pair data.frame column
  "X0",
  # plot_phyloseq example dataset
  "esophagus",
  # dcast.data.table (data.table function used in import_uparse)
  "dcast.data.table"
))
###############################################
