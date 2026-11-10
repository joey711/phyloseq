#' Plot sample abundance, richness and zero-fraction diagnostics
#'
#' Uses original IDs without melting a count matrix. Empty samples remain
#' visible; negative or non-finite abundances and ambiguous IDs are rejected.
#' This plot describes the selected nonnegative assay; it does not establish
#' that a normalization method is statistically appropriate.
#'
#' @inheritParams check_composition
#' @return A ggplot object with a panel for total abundance, observed taxa and
#'   zero fraction. Its data are available through the ggplot data element.
#' @export
plot_composition_qc <- function(x, assay="counts") {
  values <- .composition_matrix(x, assay)
  issues <- .composition_id_issues(values)
  if (length(issues)) stop(paste(issues, collapse="\n"))
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("QC plots require finite nonnegative abundances; inspect check_composition() first.")
  }
  report <- check_composition(x, assay)$samples
  if (any(!is.finite(report$Total))) {
    stop("Abundance totals overflow; rescale the data before plotting.")
  }
  metrics <- c("Total abundance", "Observed taxa", "Zero fraction")
  data <- data.frame(Sample=rep(report$Sample, 3),
                     Metric=rep(metrics, each=nrow(report)),
                     Value=c(report$Total, report$Observed, report$ZeroFraction),
                     stringsAsFactors=FALSE)
  data$Sample <- factor(data$Sample, levels=report$Sample)
  data$Metric <- factor(data$Metric, levels=metrics)
  ggplot2::ggplot(data, ggplot2::aes(x=Sample, y=Value)) +
    ggplot2::geom_col(fill="#287C8E") +
    ggplot2::facet_wrap(~Metric, scales="free_y", ncol=1) +
    ggplot2::theme_bw() + ggplot2::labs(x="Sample", y=NULL) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))
}

#' Plot PCA of CLR coordinates in Aitchison geometry
#'
#' PCA uses centered, unscaled sample CLR vectors. No additional per-feature
#' standardization or implicit zero replacement is applied. Axis percentages
#' describe variance in this CLR geometry, not statistical significance.
#'
#' @inheritParams check_composition
#' @param colour Optional sample-metadata column for point colour. Metadata
#'   are aligned by sample ID, not their current row order.
#' @return A ggplot object. Its data contain Sample, Axis1 and Axis2, and Colour
#'   when requested. Axis signs are arbitrary, as with any PCA.
#' @export
plot_composition_ordination <- function(x, assay="clr", colour=NULL) {
  check_composition(x, assay, scale="clr", stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  if (ncol(values) < 2L) stop("Ordination requires at least two samples.")
  fit <- stats::prcomp(t(values), center=TRUE, scale.=FALSE)
  variance <- fit$sdev ^ 2
  if (!is.finite(sum(variance)) || sum(variance) <= 0) {
    stop("Ordination requires variation between samples.")
  }
  data <- data.frame(Sample=colnames(values), Axis1=fit$x[, 1], Axis2=fit$x[, 2],
                     row.names=NULL, stringsAsFactors=FALSE)
  if (!is.null(colour)) {
    .composition_name(colour, "colour")
    if (inherits(x, "SummarizedExperiment")) {
      metadata <- as.data.frame(SummarizedExperiment::colData(x))
    } else if (inherits(x, "phyloseq") && !is.null(sample_data(x, FALSE))) {
      metadata <- as(sample_data(x), "data.frame")
    } else {
      stop("Sample metadata are required when `colour` is specified.")
    }
    if (sum(names(metadata) == colour) != 1L) stop("Expected one metadata column named: ", colour)
    data$Colour <- metadata[data$Sample, colour]
    mapping <- ggplot2::aes(x=Axis1, y=Axis2, colour=Colour)
  } else {
    mapping <- ggplot2::aes(x=Axis1, y=Axis2)
  }
  percentage <- 100 * variance / sum(variance)
  ggplot2::ggplot(data, mapping) + ggplot2::geom_point(size=3) + ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::labs(x=sprintf("PC1 (%.1f%%)", percentage[1]),
                  y=sprintf("PC2 (%.1f%%)", percentage[2]), colour=colour)
}

utils::globalVariables(c("Metric", "Value", "Axis1", "Axis2", "Colour"))
