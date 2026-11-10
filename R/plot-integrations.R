# Optional plotting backends share validation and ID-based annotation alignment.
.plot_require <- function(package) {
  minimum <- c(microViz="0.12.6", tidyplots="0.4.0", pheatmap="1.0.12")[[package]]
  if (!requireNamespace(package, quietly=TRUE,
                        versionCheck=list(op=">=", version=minimum))) {
    hint <- if (package == "microViz") {
      'See https://david-barnett.github.io/microViz/ for installation instructions.'
    } else {
      paste0('Install it with install.packages("', package, '").')
    }
    stop("This plot requires the optional package '", package, "' (>= ", minimum, "). ", hint)
  }
}

.plot_annotations <- function(x, margin) {
  if (inherits(x, "phyloseq")) {
    component <- if (margin == "sample") sample_data(x, FALSE) else tax_table(x, FALSE)
    if (is.null(component)) return(NULL)
    if (margin == "sample") return(as(component, "data.frame"))
    return(as.data.frame(as(component, "matrix"), stringsAsFactors=FALSE))
  }
  if (inherits(x, "SummarizedExperiment")) {
    component <- if (margin == "sample") SummarizedExperiment::colData(x) else
      SummarizedExperiment::rowData(x)
    return(as.data.frame(component))
  }
  NULL
}

.plot_fields <- function(data, fields, ids, argument) {
  if (is.null(fields)) fields <- character()
  if (!is.character(fields) || anyNA(fields) || any(!nzchar(fields)) || anyDuplicated(fields)) {
    stop("`", argument, "` must contain unique, non-empty column names.")
  }
  if (!length(fields)) return(data.frame(row.names=ids))
  if (is.null(data) || anyDuplicated(names(data)) || !all(fields %in% names(data))) {
    stop("Unknown or ambiguous annotation columns in `", argument, "`.")
  }
  if (anyDuplicated(rownames(data)) || !all(ids %in% rownames(data))) {
    stop("Annotation IDs do not match the selected data.")
  }
  data <- data[ids, fields, drop=FALSE]
  supported <- vapply(data, function(column) {
    (is.numeric(column) || is.character(column) || is.logical(column) || is.factor(column)) &&
      is.null(dim(column))
  }, logical(1))
  if (!all(supported)) {
    stop("`", argument, "` requires scalar numeric, character, logical or factor columns.")
  }
  data
}

.plot_ids <- function(ids, available, argument) {
  if (is.null(ids)) return(available)
  if (!is.character(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids) || !all(ids %in% available)) {
    stop("`", argument, "` must contain unique IDs present in the input.")
  }
  ids
}

#' Export validated assay values and annotations for plotting
#'
#' Exports values unchanged in sample-major order, with explicit character IDs.
#' There is no implicit transformation, aggregation, zero replacement or removal
#' of zero entries. The entire assay is validated before selecting taxa for
#' display, so displaying a subset of valid CLR coordinates does not recompute
#' their reference composition. Such a display subset is not a new CLR assay
#' and should not be used to compute full-composition Aitchison distances.
#'
#' @inheritParams check_composition
#' @param input_scale Expected input scale, passed to check_composition. Select
#'   clr explicitly for CLR coordinates. Recorded assay provenance is checked.
#' @param taxa Optional character IDs of taxa to display, in the desired order.
#'   NULL retains all taxa. Selection occurs after validation.
#' @param sample_vars Optional sample metadata columns to include.
#' @param tax_vars Optional taxon annotation columns to include.
#' @return A data frame with Sample, Taxon and Abundance. Requested metadata
#'   columns are aligned by ID and prefixed Sample_ or Tax_. Annotation columns
#'   must be scalar numeric, character, logical or factor vectors.
#' @export
#' @examples
#' data(esophagus)
#' head(composition_data(esophagus, taxa=taxa_names(esophagus)[1:3]))
composition_data <- function(x, assay="counts",
                              input_scale=c("abundance", "counts", "proportion", "clr"),
                              taxa=NULL, sample_vars=NULL, tax_vars=NULL) {
  input_scale <- match.arg(input_scale)
  check_composition(x, assay, scale=input_scale, stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  taxa <- .plot_ids(taxa, rownames(values), "taxa")
  values <- values[taxa, , drop=FALSE]
  samples <- .plot_fields(.plot_annotations(x, "sample"), sample_vars, colnames(values), "sample_vars")
  taxonomy <- .plot_fields(.plot_annotations(x, "taxon"), tax_vars, taxa, "tax_vars")
  data <- data.frame(Sample=rep(colnames(values), each=nrow(values)),
                     Taxon=rep(taxa, times=ncol(values)), Abundance=as.vector(values),
                     stringsAsFactors=FALSE)
  for (field in names(samples)) {
    data[[paste0("Sample_", field)]] <- samples[[field]][rep(seq_len(ncol(values)), each=nrow(values))]
  }
  for (field in names(taxonomy)) {
    data[[paste0("Tax_", field)]] <- rep(taxonomy[[field]], times=ncol(values))
  }
  data
}

#' Create an extensible tidyplots stacked abundance plot
#'
#' Builds a tidyplot with an absolute stacked bar layer on the selected values.
#' Use a precomputed proportion assay for relative compositions. No closure is
#' applied after choosing taxa: hidden taxa do not inflate displayed fractions.
#' CLR coordinates are unsuitable for compositional stacked bars and are rejected.
#'
#' @inheritParams composition_data
#' @param input_scale Expected nonnegative input scale: abundance, counts or
#'   proportion. CLR coordinates are not accepted for compositional stacked bars.
#' @param ... Plot size and style options passed to tidyplots::tidyplot. Mappings
#'   for x, y, color and fill are fixed; use composition_data for other mappings.
#' @return A tidyplot object, extendable with tidyplots functions.
#' @export
#' @examples
#' if (requireNamespace("tidyplots", quietly=TRUE)) {
#'   data(esophagus)
#'   p <- plot_composition_tidy(esophagus, taxa=taxa_names(esophagus)[1:3])
#' }
plot_composition_tidy <- function(x, assay="counts",
                                   input_scale=c("abundance", "counts", "proportion"),
                                   taxa=NULL, sample_vars=NULL, tax_vars=NULL, ...) {
  input_scale <- match.arg(input_scale)
  .plot_require("tidyplots")
  dots <- list(...)
  if ((is.null(names(dots)) && length(dots)) ||
      any(!nzchar(names(dots))) || anyDuplicated(names(dots)) ||
      any(names(dots) %in% c("data", "x", "y", "color", "colour", "fill"))) {
    stop("Use named tidyplot style options; aesthetic mappings are fixed.")
  }
  data <- composition_data(x, assay, input_scale, taxa, sample_vars, tax_vars)
  data$Sample <- factor(data$Sample, levels=unique(data$Sample))
  data$Taxon <- factor(data$Taxon, levels=unique(data$Taxon))
  plot <- tidyplots::tidyplot(data, x=Sample, y=Abundance, color=Taxon, ...)
  tidyplots::add_barstack_absolute(plot)
}

#' Draw an annotated pheatmap from a validated assay
#'
#' Validates the full assay, including recorded provenance, before optional
#' display-only taxon selection. Values are not transformed or renormalized.
#' In particular, a CLR display subset retains its original reference composition.
#' pheatmap defaults to scale="none"; an explicit scale="row" or "column" in ...
#' is a display standardization and changes the geometry used for clustering.
#'
#' @inheritParams composition_data
#' @param cluster_rows,cluster_cols NULL (default) enables clustering when at
#'   least two rows/columns are displayed. Otherwise logical values or hclust
#'   objects, as supported by pheatmap.
#' @param ... Other pheatmap::pheatmap options, including silent, filename,
#'   color, scale and clustering distances. mat, annotation_col and
#'   annotation_row are managed by this wrapper and cannot be overridden.
#' @return The native pheatmap result containing dendrograms and a gtable.
#' @export
#' @examples
#' if (requireNamespace("pheatmap", quietly=TRUE)) {
#'   data(esophagus)
#'   h <- plot_composition_heatmap(esophagus, taxa=taxa_names(esophagus)[1:5],
#'                                  silent=TRUE)
#' }
plot_composition_heatmap <- function(x, assay="counts",
                                      input_scale=c("abundance", "counts", "proportion", "clr"),
                                      taxa=NULL, sample_vars=NULL, tax_vars=NULL,
                                      cluster_rows=NULL, cluster_cols=NULL, ...) {
  input_scale <- match.arg(input_scale)
  .plot_require("pheatmap")
  dots <- list(...)
  if ((is.null(names(dots)) && length(dots)) || any(!nzchar(names(dots))) ||
      anyDuplicated(names(dots)) || any(names(dots) %in% c("mat", "annotation_col", "annotation_row", "annotation"))) {
    stop("Use named pheatmap options; values and annotations are managed by the wrapper.")
  }
  check_composition(x, assay, scale=input_scale, stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  taxa <- .plot_ids(taxa, rownames(values), "taxa")
  values <- values[taxa, , drop=FALSE]
  samples <- .plot_fields(.plot_annotations(x, "sample"), sample_vars, colnames(values), "sample_vars")
  taxonomy <- .plot_fields(.plot_annotations(x, "taxon"), tax_vars, taxa, "tax_vars")
  # pheatmap accepts discrete factors/strings or continuous numeric annotations.
  for (field in names(samples)) if (is.logical(samples[[field]])) samples[[field]] <- factor(samples[[field]])
  for (field in names(taxonomy)) if (is.logical(taxonomy[[field]])) taxonomy[[field]] <- factor(taxonomy[[field]])
  if (is.null(cluster_rows)) cluster_rows <- nrow(values) > 1L
  if (is.null(cluster_cols)) cluster_cols <- ncol(values) > 1L
  # Constant assays need nondegenerate breaks, including a one-cell heatmap.
  if (all(values == values[1L]) && is.null(dots$breaks)) {
    limits <- range(c(0, values[1L]))
    if (limits[1L] == limits[2L]) limits <- c(-1, 1)
    ncolors <- if (is.null(dots$color)) 100L else length(dots$color)
    dots$breaks <- seq(limits[1L], limits[2L], length.out=ncolors + 1L)
  }
  do.call(pheatmap::pheatmap,
          c(list(mat=values, annotation_col=if (ncol(samples)) samples else NA,
                 annotation_row=if (ncol(taxonomy)) taxonomy else NA,
                 cluster_rows=cluster_rows, cluster_cols=cluster_cols), dots))
}

#' Adapt a validated assay for the microViz workflow
#'
#' Creates a microViz psExtra with tax_transform(trans="identity", rank="unique",
#' keep_counts=FALSE). No additional numerical transformation or zero replacement
#' occurs. The identity label means that this adapter leaves values unchanged,
#' not that a selected CLR/proportion assay contains raw counts. Use explicit
#' ordination methods and appropriate downstream plotting options for that scale.
#'
#' phyloseq inputs retain their components. Experiment inputs supply the selected
#' assay, scalar sample metadata and character/factor row annotations as taxonomy.
#' Other assays, experiment metadata, row trees and reference sequences are not
#' transferred to the plotting adapter. They remain in the unchanged source
#' experiment, which is the authoritative object for provenance validation.
#' A missing taxonomy table is supplied with feature IDs. microViz adds a unique
#' taxonomic rank for feature IDs. Its subsequent operations
#' do not participate in phyloseqV2's assay fingerprint checks.
#'
#' @inheritParams composition_data
#' @return A microViz psExtra, usable with microViz::comp_barplot,
#'   microViz::comp_heatmap, microViz::ord_calc and other microViz functions.
#' @export
#' @examples
#' if (requireNamespace("microViz", quietly=TRUE)) {
#'   data(esophagus)
#'   mv <- as_microviz(otu_table(esophagus))
#'   p <- microViz::comp_barplot(mv, tax_level="unique", sample_order="default")
#' }
as_microviz <- function(x, assay="counts",
                         input_scale=c("abundance", "counts", "proportion", "clr")) {
  input_scale <- match.arg(input_scale)
  .plot_require("microViz")
  check_composition(x, assay, scale=input_scale, stop_on_error=TRUE)
  values <- .composition_matrix(x, assay)
  if (inherits(x, "phyloseq")) {
    ps <- x
  } else {
    components <- list(otu_table(values, taxa_are_rows=TRUE))
    samples <- .plot_annotations(x, "sample")
    samples <- .plot_fields(samples, names(samples), colnames(values), "sample metadata")
    if (ncol(samples)) components <- c(components, list(sample_data(samples)))
    taxonomy <- .plot_annotations(x, "taxon")
    fields <- names(taxonomy)[vapply(taxonomy, function(column) is.character(column) || is.factor(column), logical(1))]
    taxonomy <- .plot_fields(taxonomy, fields, rownames(values), "taxon annotations")
    if (!ncol(taxonomy)) {
      taxonomy <- data.frame(unique=rownames(values), row.names=rownames(values),
                              stringsAsFactors=FALSE)
    }
    components <- c(components, list(tax_table(as.matrix(taxonomy))))
    ps <- do.call(phyloseq, components)
  }
  if (is.null(tax_table(ps, FALSE))) {
    tax_table(ps) <- tax_table(matrix(taxa_names(ps), ncol=1,
                                      dimnames=list(taxa_names(ps), "unique")))
  }
  result <- microViz::tax_transform(ps, trans="identity", rank="unique", keep_counts=FALSE)
  adapted <- .composition_matrix(microViz::ps_get(result))
  if (!identical(dimnames(adapted), dimnames(values)) || !all(adapted == values)) {
    stop("microViz adaptation changed values or identifiers; refusing altered data.")
  }
  result
}

utils::globalVariables(c("Sample", "Taxon", "Abundance"))
