.silva_version <- function(version) {
  if (!is.character(version) || length(version) != 1L || is.na(version) ||
      !grepl("^[0-9]+([.][0-9]+)?$", version)) stop("`version` must be an explicit SILVA release string, e.g. \"144\".")
}

.silva_check <- function(reference) {
  if (!inherits(reference, "silva_reference")) stop("Expected a silva_reference descriptor.")
  for (field in c("fasta", "quality")) {
    file <- reference[[field]]
    if (!is.character(file) || length(file) != 1L || !file.exists(file)) stop("Missing reference file: ", field)
    if (!identical(unname(tools::md5sum(file)), reference$checksums[[field]])) {
      stop("Reference file changed or failed its checksum: ", field)
    }
  }
}

.silva_md5 <- function(file) {
  lines <- readLines(file, warn=FALSE)
  words <- unlist(strsplit(lines, "[[:space:]]+"), use.names=FALSE)
  hashes <- words[grepl("^[0-9a-fA-F]{32}$", words)]
  if (length(hashes) != 1L) stop("Expected exactly one MD5 checksum in: ", file)
  tolower(hashes)
}

#' Obtain a versioned SILVA SSU Ref NR99 reference
#'
#' Uses explicit releases and the official unaligned, gene-trimmed FASTA plus
#' quality exports. SSU contains both 16S and 18S: use filter_silva_reference to
#' select bacterial/archaeal 16S. A gene-trimmed export is not an amplicon-only
#' database and is not proof that every sequence is complete.
#'
#' @param version Explicit release string. Official downloads currently support
#'   "138.2" and "144"; other releases can be supplied as local files.
#' @param cache_dir Writable cache directory; defaults to a directory under
#'   the current R temporary directory. Choose a persistent directory to reuse it.
#' @param download If TRUE, download uncached official files and verify publisher
#'   MD5 sidecars. Downloads can be hundreds of MB. FALSE never accesses the network.
#' @param fasta,quality Optional paired local FASTA and SILVA quality-export paths.
#'   Local files receive computed fingerprints, not publisher-verification claims.
#' @return A silva_reference descriptor with paths, source, release, checksums
#'   and publisher_verified. Sequence data are not loaded into memory.
#' @export
silva_reference <- function(version, cache_dir=file.path(tempdir(), "phyloseqV2-silva"),
                             download=FALSE, fasta=NULL, quality=NULL) {
  .silva_version(version)
  if (!is.logical(download) || length(download) != 1L || is.na(download)) stop("`download` must be TRUE or FALSE.")
  if (download) {
    # The base 60-second timeout can be too short for a full reference database.
    old_options <- options(timeout=max(600, getOption("timeout", 60)))
    on.exit(options(old_options), add=TRUE)
  }
  local <- !is.null(fasta) || !is.null(quality)
  if (local) {
    if (is.null(fasta) || is.null(quality) || download) stop("Supply both local files and use download=FALSE.")
    for (file in list(fasta, quality)) {
      .composition_name(file, "local file")
      if (!file.exists(file) || dir.exists(file)) stop("Local reference file does not exist: ", file)
    }
    source <- "local"
    verified <- FALSE
    urls <- NULL
  } else {
    if (!version %in% c("138.2", "144")) stop("Official download supports 138.2 and 144; supply local files for other releases.")
    .composition_name(cache_dir, "cache_dir")
    directory <- file.path(path.expand(cache_dir), paste0("SILVA_", version, "_SSURef_NR99"))
    filenames <- c(fasta=paste0("SILVA_", version, "_SSURef_NR99_tax_silva_trunc.fasta.gz"),
                   quality=paste0("SILVA_", version, "_SSURef_Nr99.quality.gz"))
    base <- paste0("https://www.arb-silva.de/fileadmin/silva_databases/release_",
                   gsub(".", "_", version, fixed=TRUE), "/Exports/")
    urls <- c(fasta=paste0(base, filenames[["fasta"]]),
               quality=paste0(base, "quality/", filenames[["quality"]]))
    paths <- file.path(directory, filenames)
    names(paths) <- names(filenames)
    staging_directories <- character()
    on.exit(unlink(staging_directories, recursive=TRUE), add=TRUE)
    for (field in names(paths)) {
      destination <- paths[[field]]
      sidecar <- paste0(destination, ".md5")
      if (!file.exists(destination) || !file.exists(sidecar)) {
        if (!download) stop("SILVA reference is not cached. Use download=TRUE or supply local fasta/quality files.")
        if (!dir.exists(directory) && !dir.create(directory, recursive=TRUE)) stop("Cannot create SILVA cache directory.")
        stage <- tempfile("download-", tmpdir=directory)
        dir.create(stage)
        staging_directories <- c(staging_directories, stage)
        staged_file <- file.path(stage, basename(destination))
        utils::download.file(paste0(urls[[field]], ".md5"), paste0(staged_file, ".md5"), mode="wb", quiet=TRUE)
        expected <- .silva_md5(paste0(staged_file, ".md5"))
        if (file.exists(destination)) {
          if (!identical(unname(tools::md5sum(destination)), expected)) stop("Existing cache file has an invalid checksum; use a fresh cache directory.")
        } else {
          utils::download.file(urls[[field]], staged_file, mode="wb", quiet=TRUE)
          if (!identical(unname(tools::md5sum(staged_file)), expected)) stop("Downloaded SILVA file failed its publisher checksum.")
          if (!file.rename(staged_file, destination)) stop("Cannot publish the verified SILVA cache file.")
        }
        if (!file.rename(paste0(staged_file, ".md5"), sidecar)) stop("Cannot publish the checksum sidecar.")
      }
      if (!identical(unname(tools::md5sum(destination)), .silva_md5(sidecar))) stop("Cached SILVA file failed its publisher checksum.")
    }
    if (download) {
      for (file in c("LICENSE.txt", "CITATION.txt")) {
        destination <- file.path(directory, file)
        if (!file.exists(destination)) {
          url <- if (file == "LICENSE.txt") paste0(base, file) else paste0(sub("Exports/$", "", base), file)
          utils::download.file(url, destination, mode="wb", quiet=TRUE)
        }
      }
    }
    fasta <- paths[["fasta"]]
    quality <- paths[["quality"]]
    source <- "official-cache"
    verified <- TRUE
  }
  result <- list(version=version, dataset="SSURef_NR99", target="16S", source=source,
                  publisher_verified=verified, urls=urls,
                  fasta=normalizePath(fasta), quality=normalizePath(quality),
                  checksums=list(fasta=unname(tools::md5sum(fasta)), quality=unname(tools::md5sum(quality))))
  structure(result, class="silva_reference")
}

.silva_quality <- function(file) {
  # Read identifiers as characters: accession and gene coordinates jointly key a
  # reference. Do not join on accession alone or type-convert leading-zero IDs.
  connection <- gzfile(file, "rt")
  header <- readLines(connection, n=1L, warn=FALSE)
  close(connection)
  if (!length(header)) stop("The SILVA quality file is empty.")
  columns <- strsplit(header, "\t", fixed=TRUE)[[1L]]
  required <- c("Primary Accession", "Start", "Stop")
  if (anyDuplicated(columns) || !all(required %in% columns)) stop("Unsupported SILVA quality schema: missing or duplicate identifier columns.")
  data <- utils::read.delim(file, check.names=FALSE, quote="", comment.char="",
                            colClasses=c("Primary Accession"="character", Start="character", Stop="character"),
                            na.strings=c("", "NA"), stringsAsFactors=FALSE)
  if (anyDuplicated(names(data)) || !all(required %in% names(data))) stop("Unsupported SILVA quality schema: missing or duplicate identifier columns.")
  if (anyNA(data[, required]) || any(!nzchar(data[[required[1L]]])) ||
      any(!grepl("^[0-9]+$", data$Start)) || any(!grepl("^[0-9]+$", data$Stop))) stop("Invalid SILVA quality identifiers.")
  ids <- paste(data[[required[1L]]], data$Start, data$Stop, sep=".")
  .asv_ids(ids, "SILVA reference IDs in quality metadata")
  rownames(data) <- ids
  data
}

.silva_chunks <- function(file, chunk_size, callback) {
  input <- XVector::open_input_files(file)
  repeat {
    chunk <- Biostrings::readBStringSet(input, nrec=chunk_size)
    if (!length(chunk)) break
    callback(chunk)
  }
  invisible(NULL)
}

#' Filter SILVA references for bacterial and archaeal 16S phylogeny
#'
#' Streams the unaligned FASTA and writes retained sequences plus a gzipped
#' per-reference QC report. SILVA quality metadata are loaded in memory, but the
#' complete sequence database/alignment is not. Default length minima follow the
#' SILVA Ref domain distinction; they do not prove complete genes. Quality cutoffs
#' are explicit, editable defaults, not universally optimal scientific choices.
#'
#' Scores are on a 0-100 scale; ambiguity, homopolymer and vector thresholds are
#' fractions. Computed ambiguity and SILVA-reported ambiguity must both pass.
#' Metadata are keyed on accession.start.stop. Requested absent fields are errors,
#' including Pintail in release 144. No score or Phred quality is invented.
#'
#' @param reference A silva_reference descriptor from silva_reference.
#' @param output_dir New/nonconflicting output directory. Existing result files
#'   are never overwritten. Raw reference files remain unchanged.
#' @param domains Subset of c("Bacteria", "Archaea"). Eukaryotic 18S is excluded.
#' @param min_length Named minimum ungapped lengths for Bacteria and Archaea.
#' @param max_length Maximum ungapped sequence length.
#' @param max_ambiguous,max_homopolymer,max_vector Maximum fractions, from 0 to 1.
#' @param min_sequence_quality,min_alignment_quality Minimum SILVA scores, 0-100.
#' @param min_pintail Optional minimum historical Pintail score. NULL does not
#'   apply that filter and does not claim that a chimera check was performed.
#' @param missing_quality "exclude" (default) reports/excludes missing values;
#'   "error" stops when a retained-domain sequence lacks required quality values.
#' @param chunk_size Number of FASTA records read at once.
#' @return A silva_filtered_reference descriptor with filtered FASTA, source
#'   fingerprints, thresholds, retained/excluded counts and qc_report path.
#'   A result with zero retained references is valid for inspecting exclusions;
#'   it cannot be used for phylogenetic inference.
#' @export
filter_silva_reference <- function(reference, output_dir, domains=c("Bacteria", "Archaea"),
                                    min_length=c(Bacteria=1200, Archaea=900), max_length=3000,
                                    max_ambiguous=0.02, max_homopolymer=0.02, max_vector=0.02,
                                    min_sequence_quality=95, min_alignment_quality=75, min_pintail=NULL,
                                    missing_quality=c("exclude", "error"), chunk_size=10000L) {
  .silva_check(reference)
  missing_quality <- match.arg(missing_quality)
  if (!is.character(domains) || !length(domains) || anyNA(domains) || anyDuplicated(domains) ||
      !all(domains %in% c("Bacteria", "Archaea"))) stop("`domains` must select Bacteria and/or Archaea.")
  if (!is.numeric(min_length) || is.complex(min_length) || anyNA(min_length) ||
      any(!is.finite(min_length)) || any(min_length < 1) || any(min_length != floor(min_length)) ||
      anyDuplicated(names(min_length)) || !all(domains %in% names(min_length))) stop("`min_length` must provide named, positive whole lengths for selected domains.")
  .asv_scalar(max_length, "max_length", max(min_length[domains]), whole=TRUE)
  for (argument in c("max_ambiguous", "max_homopolymer", "max_vector")) .asv_scalar(get(argument), argument, 0, 1)
  for (argument in c("min_sequence_quality", "min_alignment_quality")) .asv_scalar(get(argument), argument, 0, 100)
  if (!is.null(min_pintail)) .asv_scalar(min_pintail, "min_pintail", 0, 100)
  if (reference$version == "144" && !is.null(min_pintail)) stop("Pintail was discontinued in SILVA 144; choose available quality criteria.")
  .asv_scalar(chunk_size, "chunk_size", 1, whole=TRUE)
  quality <- .silva_quality(reference$quality)
  fields <- c(SequenceQuality="Sequence Quality", AlignmentQuality="Alignment Quality",
               ReportedAmbiguousFraction="% Ambiguities", HomopolymerFraction="% Homopolymers",
               VectorFraction="% Vector Contamination")
  if (!is.null(min_pintail)) fields <- c(fields, Pintail="Pintail Quality")
  if (!all(fields %in% names(quality))) stop("Required quality fields are unavailable for this release: ",
                                           paste(setdiff(fields, names(quality)), collapse=", "))
  for (field in fields) {
    before <- quality[[field]]
    values <- suppressWarnings(as.numeric(before))
    if (any(!is.na(before) & (is.na(values) | !is.finite(values) | values < 0 | values > 100))) {
      stop("Invalid 0-100 quality values in field: ", field)
    }
    quality[[field]] <- values
  }
  .composition_name(output_dir, "output_dir")
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  output <- file.path(output_dir, c(fasta="reference-16S.fasta", report="reference-qc.tsv.gz", manifest="reference.rds"))
  names(output) <- c("fasta", "report", "manifest")
  if (any(file.exists(output))) stop("Reference output files already exist; choose a new output directory.")
  stage <- tempfile("silva-qc-", tmpdir=output_dir)
  dir.create(stage)
  fasta <- file.path(stage, "reference.fasta")
  file.create(fasta)
  report_file <- file.path(stage, "report.tsv.gz")
  report_connection <- gzfile(report_file, "wt")
  report_open <- TRUE
  on.exit({
    if (report_open) close(report_connection)
    unlink(stage, recursive=TRUE)
  }, add=TRUE)
  seen <- logical(nrow(quality))
  unknown <- new.env(hash=TRUE, parent=emptyenv())
  total <- retained <- 0L
  first <- TRUE
  reason_counts <- integer()
  .silva_chunks(reference$fasta, chunk_size, function(chunk) {
    headers <- names(chunk)
    ids <- sub("[[:space:]].*$", "", headers)
    .asv_ids(ids, "SILVA FASTA IDs within a chunk")
    taxonomy <- ifelse(grepl("[[:space:]]", headers), sub("^[^[:space:]]+[[:space:]]+", "", headers), "")
    domain <- sub(";.*$", "", taxonomy)
    positions <- match(ids, rownames(quality))
    known <- !is.na(positions)
    if (any(seen[positions[known]])) stop("Duplicate reference IDs across FASTA chunks.")
    seen[positions[known]] <<- TRUE
    for (id in ids[!known]) {
      if (exists(id, envir=unknown, inherits=FALSE)) stop("Duplicate reference ID: ", id)
      assign(id, TRUE, envir=unknown)
    }
    names(chunk) <- ids
    dna <- .asv_dna(chunk, unambiguous=FALSE)
    lengths <- Biostrings::width(dna)
    ambiguous <- (lengths - rowSums(Biostrings::letterFrequency(dna, c("A", "C", "G", "T")))) / lengths
    report <- data.frame(ReferenceID=ids, Domain=domain, Length=lengths,
                           AmbiguousFraction=ambiguous, stringsAsFactors=FALSE)
    for (field in names(fields)) {
      report[[field]] <- quality[[fields[[field]]]][positions]
      if (grepl("Fraction$", field)) report[[field]] <- report[[field]] / 100
    }
    missing <- !known | !stats::complete.cases(report[, names(fields), drop=FALSE])
    if (missing_quality == "error" && any(missing & domain %in% domains)) stop("Required quality values are missing for retained-domain reference sequences.")
    failures <- list(domain=!domain %in% domains,
                      length=lengths < unname(min_length[domain]) | lengths > max_length,
                      missing_quality=missing,
                      ambiguity=ambiguous > max_ambiguous | report$ReportedAmbiguousFraction > max_ambiguous,
                      sequence_quality=report$SequenceQuality < min_sequence_quality,
                      alignment_quality=report$AlignmentQuality < min_alignment_quality,
                      homopolymer=report$HomopolymerFraction > max_homopolymer,
                      vector=report$VectorFraction > max_vector)
    if (!is.null(min_pintail)) failures$pintail <- report$Pintail < min_pintail
    reasons <- rep("", length(dna))
    for (reason in names(failures)) {
      failing <- failures[[reason]]
      failing[is.na(failing)] <- FALSE
      reasons[failing] <- ifelse(nzchar(reasons[failing]), paste(reasons[failing], reason, sep=";"), reason)
      previous <- reason_counts[reason]
      if (is.na(previous)) previous <- 0L
      reason_counts[reason] <<- previous + sum(failing)
    }
    keep <- !nzchar(reasons)
    report$Keep <- keep
    report$Reason <- reasons
    names(dna) <- headers
    if (any(keep)) Biostrings::writeXStringSet(dna[keep], fasta, append=TRUE)
    utils::write.table(report, report_connection, sep="\t", quote=FALSE, row.names=FALSE, col.names=first, na="NA")
    first <<- FALSE
    total <<- total + length(dna)
    retained <<- retained + sum(keep)
  })
  if (!total) stop("The input reference FASTA contains no sequences.")
  close(report_connection)
  report_open <- FALSE
  parameters <- list(domains=domains, min_length=min_length, max_length=max_length,
                       max_ambiguous=max_ambiguous, max_homopolymer=max_homopolymer, max_vector=max_vector,
                       min_sequence_quality=min_sequence_quality, min_alignment_quality=min_alignment_quality,
                       min_pintail=min_pintail, missing_quality=missing_quality)
  if (!file.rename(fasta, output[["fasta"]]) || !file.rename(report_file, output[["report"]])) stop("Cannot publish filtered reference files.")
  result <- reference
  result$source_reference <- reference
  result$source <- "filtered"
  result$publisher_verified <- FALSE
  result$fasta <- normalizePath(output[["fasta"]])
  result$checksums$fasta <- unname(tools::md5sum(result$fasta))
  result$qc_report <- normalizePath(output[["report"]])
  result$qc_report_md5 <- unname(tools::md5sum(result$qc_report))
  result$parameters <- parameters
  result$summary <- list(total=total, retained=retained, excluded=total-retained, failure_counts=reason_counts)
  class(result) <- c("silva_filtered_reference", "silva_reference")
  saveRDS(result, output[["manifest"]])
  result
}
