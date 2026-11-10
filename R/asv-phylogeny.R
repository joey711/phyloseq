.asv_reference_chunks <- function(reference, callback, chunk_size=10000L) {
  if (inherits(reference, "silva_reference")) {
    .silva_check(reference)
    if (!inherits(reference, "silva_filtered_reference")) stop("Filter the SILVA reference before matching or inference.")
    .silva_chunks(reference$fasta, chunk_size, function(chunk) {
      names(chunk) <- sub("[[:space:]].*$", "", names(chunk))
      callback(.asv_dna(chunk, unambiguous=FALSE))
    })
  } else {
    callback(.asv_dna(reference, unambiguous=FALSE))
  }
  invisible(NULL)
}

#' Match observed ASVs to quality-filtered reference genes
#'
#' Performs a full-query, ungapped subsequence search using Biostrings. This
#' conservative candidate search is exact by default, with an explicit optional
#' mismatch allowance; it is not a general gapped homology search or a taxonomic
#' classifier. Identity is matches divided by query length, and query coverage
#' is 1. Reference coverage is not required to be 1 for a short ASV.
#'
#' Retains one best position per reference and strand, preserving one-to-many
#' reference matches. Reports no-hit queries and candidate truncation. No observed
#' sequence is replaced by a full-length reference or silently removed.
#'
#' @param x An input accepted by asv_sequences.
#' @param reference A filtered SILVA reference or a named XStringSet of already
#'   curated homologous reference sequences. Manually supplied sequences bypass
#'   SILVA QC and are labelled as caller-curated.
#' @param max_mismatches Maximum mismatches across the entire query. No indels.
#' @param both_strands Search both observed sequence and its reverse complement.
#' @param max_hits Maximum returned candidates per ASV, or Inf to retain all.
#'   Candidates are ranked by mismatches, reference ID, strand and position;
#'   truncated searches are marked, not reported as unique placements.
#' @param chunk_size FASTA records per chunk for on-disk references.
#' @return An asv_reference_matches list with query_sequences, reference,
#'   matches, per-query summary, parameters and query fingerprints. Matching
#'   scores are not placement probabilities or evidence of species identity.
#' @export
#' @examples
#' queries <- c(asv1="ACGTAC")
#' refs <- c(ref1="TTACGTACTT", ref2="GGACGTACAA")
#' match_asvs_to_reference(queries, refs)$matches
match_asvs_to_reference <- function(x, reference, max_mismatches=0L, both_strands=TRUE,
                                     max_hits=100L, chunk_size=10000L) {
  queries <- asv_sequences(x)
  .asv_scalar(max_mismatches, "max_mismatches", 0, min(Biostrings::width(queries)) - 1, whole=TRUE)
  .asv_scalar(chunk_size, "chunk_size", 1, whole=TRUE)
  if (!identical(max_hits, Inf)) .asv_scalar(max_hits, "max_hits", 1, whole=TRUE)
  if (!is.logical(both_strands) || length(both_strands) != 1L || is.na(both_strands)) stop("`both_strands` must be TRUE or FALSE.")
  empty <- data.frame(QueryID=character(), ReferenceID=character(), Strand=character(),
                       Start=integer(), End=integer(), Mismatches=integer(), Identity=numeric(),
                       QueryCoverage=numeric(), stringsAsFactors=FALSE)
  candidates <- stats::setNames(rep(list(empty), length(queries)), names(queries))
  totals <- integer(length(queries))
  best <- rep(Inf, length(queries))
  orientation <- rep(0L, length(queries))
  .asv_reference_chunks(reference, function(subjects) {
    # Filtered descriptors were checked for global ID uniqueness during QC;
    # in-memory references are checked by .asv_dna. Avoid a database-sized hash.
    for (query_index in seq_along(queries)) {
      strands <- if (both_strands) c("+", "-") else "+"
      found <- list()
      for (strand in strands) {
        query <- queries[[query_index]]
        if (strand == "-") query <- Biostrings::reverseComplement(query)
        locations <- Biostrings::vmatchPattern(query, subjects, max.mismatch=max_mismatches, fixed=TRUE)
        reference_indices <- which(S4Vectors::elementNROWS(locations) > 0L)
        for (index in reference_indices) {
          starts <- BiocGenerics::start(locations[[index]])
          mismatch <- vapply(starts, function(start) {
            window <- Biostrings::subseq(subjects[[index]], start=start, width=length(query))
            sum(as.integer(window) != as.integer(query))
          }, integer(1))
          chosen <- which.min(mismatch)
          n <- mismatch[chosen]
          found[[length(found) + 1L]] <- data.frame(
            QueryID=names(queries)[query_index], ReferenceID=names(subjects)[index],
            Strand=strand, Start=starts[chosen], End=starts[chosen] + length(query) - 1L,
            Mismatches=n, Identity=1 - n / length(query), QueryCoverage=1,
            stringsAsFactors=FALSE)
        }
      }
      if (!length(found)) next
      data <- do.call(rbind, found)
      totals[query_index] <<- totals[query_index] + nrow(data)
      current_best <- min(data$Mismatches)
      bits <- 0L
      if (any(data$Strand[data$Mismatches == current_best] == "+")) bits <- bitwOr(bits, 1L)
      if (any(data$Strand[data$Mismatches == current_best] == "-")) bits <- bitwOr(bits, 2L)
      if (current_best < best[query_index]) orientation[query_index] <<- bits
      if (current_best == best[query_index]) orientation[query_index] <<- bitwOr(orientation[query_index], bits)
      best[query_index] <<- min(best[query_index], current_best)
      all <- rbind(candidates[[query_index]], data)
      all <- all[order(all$Mismatches, all$ReferenceID, all$Strand, all$Start), , drop=FALSE]
      if (nrow(all) > max_hits) all <- all[seq_len(max_hits), , drop=FALSE]
      rownames(all) <- NULL
      candidates[[query_index]] <<- all
    }
  }, chunk_size)
  matches <- do.call(rbind, unname(candidates))
  rownames(matches) <- NULL
  summary <- data.frame(QueryID=names(queries), Candidates=totals,
                        Returned=vapply(candidates, nrow, integer(1)),
                        Truncated=totals > max_hits, Status=ifelse(totals > 0, "matched", "no_hit"),
                        BestMismatches=ifelse(is.finite(best), best, NA_real_),
                        Strand=c(NA_character_, "+", "-", "ambiguous")[orientation + 1L],
                        row.names=NULL, stringsAsFactors=FALSE)
  hashes <- vapply(as.character(queries), digest::digest, character(1), algo="xxhash64")
  structure(list(query_sequences=queries, reference=reference, matches=matches, summary=summary,
                   parameters=list(method="ungapped_full_query", max_mismatches=max_mismatches,
                                     both_strands=both_strands, max_hits=max_hits),
                   query_hashes=hashes,
                   candidate_hash=digest::digest(list(matches=matches, summary=summary), algo="xxhash64"),
                   reference_hash=if (!inherits(reference,"silva_reference")) digest::digest(reference,algo="xxhash64") else NULL),
            class="asv_reference_matches")
}

#' Infer a joint maximum-likelihood tree of ASV fragments and reference genes
#'
#' A targeted, reference-enriched joint alignment/tree, not a phylogenetic
#' placement algorithm and not reconstruction of unobserved ASV bases. DECIPHER
#' aligns the actual reference and ASV sequences; phangorn fits a joint ML tree
#' starting from neighbor joining. The resulting tree is unrooted. Do not use it
#' for UniFrac without an explicit, scientifically justified rooting decision.
#'
#' All queries must have candidate references. Ambiguous strand assignments,
#' truncated candidate sets and overly large trees require explicit resolution.
#' Full-query mismatches may miss true homologs containing indels; such no-hit
#' cases require a gapped candidate search rather than relaxed inference checks.
#' Duplicate sequences with different IDs remain separate observed features.
#'
#' @param matches Output of match_asvs_to_reference.
#' @param reference_ids Optional explicit reference-ID subset from the returned
#'   candidate table. At least three references and one hit per query are required.
#' @param query_strands Optional named vector of "+"/"-" choices for all queries.
#'   Required when candidate evidence does not resolve orientation.
#' @param model Substitution model, "GTR" or "JC". Model settings are recorded.
#' @param gamma_categories Number of discrete gamma categories; 1 disables gamma.
#' @param bootstrap Number of nonparametric ML bootstrap replicates, default 0.
#'   Zero means no support assessment, and is explicitly recorded as such.
#' @param min_overlap Minimum unambiguous alignment sites shared by every pair.
#'   Disjoint amplicon regions cannot be combined by this joint-tree method.
#' @param max_sequences Maximum total sequences for this targeted in-R backend.
#' @param seed Explicit RNG seed. The caller's RNG state is restored.
#' @param threads DECIPHER worker count. Bootstrap runs sequentially.
#' @return An asv_phylogeny list with tree (ASV/reference context), asv_tree
#'   (only query IDs, or NULL for one ASV), alignment, tip_data, fitted model,
#'   bootstrap trees, provenance and the unchanged observed query sequences.
#'   No input count table or existing container tree is overwritten.
#' @export
infer_asv_tree <- function(matches, reference_ids=NULL, query_strands=NULL,
                            model=c("GTR", "JC"), gamma_categories=4L, bootstrap=0L,
                            min_overlap=20L, max_sequences=300L, seed=1L, threads=1L) {
  if (!inherits(matches, "asv_reference_matches")) stop("Expected output from match_asvs_to_reference.")
  had_rng <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_rng) saved_rng <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had_rng) assign(".Random.seed", saved_rng, envir=.GlobalEnv) else
    if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv), add=TRUE)
  for (package in c("DECIPHER", "phangorn")) .composition_require(package)
  model <- match.arg(model)
  for (argument in c("gamma_categories", "min_overlap", "max_sequences", "threads")) .asv_scalar(get(argument), argument, 1, whole=TRUE)
  .asv_scalar(bootstrap, "bootstrap", 0, whole=TRUE)
  .asv_scalar(seed, "seed", 0, .Machine$integer.max, whole=TRUE)
  queries <- .asv_dna(matches$query_sequences)
  current_hashes <- vapply(as.character(queries), digest::digest, character(1), algo="xxhash64")
  if (!identical(current_hashes, matches$query_hashes)) stop("Observed query sequences changed after candidate matching.")
  if (!identical(digest::digest(list(matches=matches$matches,summary=matches$summary),algo="xxhash64"),matches$candidate_hash)) {
    stop("Candidate matches or their summary changed; rerun candidate matching.")
  }
  if (!inherits(matches$reference,"silva_reference") &&
      !identical(digest::digest(matches$reference,algo="xxhash64"),matches$reference_hash)) stop("Reference sequences changed after candidate matching.")
  summary <- matches$summary[match(names(queries), matches$summary$QueryID), , drop=FALSE]
  if (anyNA(summary$Status) || any(summary$Status != "matched")) stop("Some ASVs have no reference candidates; inspect the match summary.")
  if (is.null(reference_ids)) {
    if (any(summary$Truncated)) stop("Candidate sets were truncated; increase max_hits or explicitly select reference_ids.")
    reference_ids <- unique(matches$matches$ReferenceID)
  }
  .asv_ids(reference_ids, "selected reference IDs")
  if (!all(reference_ids %in% matches$matches$ReferenceID)) stop("Selected reference IDs must occur in the candidate table.")
  if (length(reference_ids) < 3L) stop("Joint inference requires at least three reference sequences.")
  represented <- unique(matches$matches$QueryID[matches$matches$ReferenceID %in% reference_ids])
  if (!all(names(queries) %in% represented)) stop("Selected references must represent every query.")
  if (length(reference_ids) + length(queries) > max_sequences) stop("Too many sequences for the targeted in-R backend; select a smaller reference context or raise max_sequences explicitly.")
  if (is.null(query_strands)) query_strands <- stats::setNames(summary$Strand, summary$QueryID)
  if (!is.character(query_strands) || anyNA(query_strands) || anyDuplicated(names(query_strands)) ||
      !setequal(names(query_strands), names(queries)) || !all(query_strands %in% c("+", "-"))) {
    stop("Resolve query orientation using named '+'/'-' query_strands for every ASV.")
  }
  compatible <- matches$matches$ReferenceID %in% reference_ids &
    matches$matches$Strand == unname(query_strands[matches$matches$QueryID])
  if (!all(names(queries) %in% matches$matches$QueryID[compatible])) {
    stop("Selected references do not support the chosen query orientation for every ASV.")
  }
  oriented <- queries
  reverse <- query_strands[names(queries)] == "-"
  oriented[reverse] <- Biostrings::reverseComplement(oriented[reverse])
  reference_parts <- list()
  .asv_reference_chunks(matches$reference, function(chunk) {
    keep <- names(chunk) %in% reference_ids
    if (any(keep)) reference_parts[[length(reference_parts) + 1L]] <<- chunk[keep]
  })
  if (!length(reference_parts)) stop("Selected references are missing from the source.")
  references <- do.call(c, reference_parts)
  if (!setequal(names(references), reference_ids) || anyDuplicated(names(references))) stop("Selected reference IDs do not match the source.")
  references <- references[reference_ids]
  internal <- c(paste0("REF", seq_along(references)), paste0("ASV", seq_along(queries)))
  tip_data <- data.frame(InternalID=internal, Origin=c(rep("reference", length(references)), rep("ASV", length(queries))),
                          FeatureID=c(names(references), names(queries)),
                          TipLabel=c(paste0("REF::", names(references)), paste0("ASV::", names(queries))),
                          stringsAsFactors=FALSE)
  combined <- c(references, oriented)
  names(combined) <- internal
  set.seed(seed)
  alignment <- DECIPHER::AlignSeqs(combined, processors=threads, verbose=FALSE)
  if (!setequal(names(alignment), internal) || anyDuplicated(names(alignment))) stop("Alignment changed sequence identifiers.")
  alignment <- alignment[internal]
  if (!identical(unname(gsub("[-.]", "", as.character(alignment))), unname(as.character(combined)))) {
    stop("Alignment changed observed nucleotide content.")
  }
  letters <- as.matrix(alignment)
  observed <- matrix(letters %in% c("A", "C", "G", "T"), nrow=nrow(letters))
  shared <- tcrossprod(observed)
  if (any(shared[lower.tri(shared)] < min_overlap)) stop("Insufficient shared alignment sites; restrict the target region or use a placement workflow.")
  data <- phangorn::phyDat(letters, type="DNA")
  distances <- phangorn::dist.ml(data)
  if (any(!is.finite(distances))) stop("Pairwise starting distances are undefined; inspect alignment overlap and saturation.")
  initial <- phangorn::NJ(distances)
  # NJ can yield negative starting edges; these are numerical initialization,
  # not the reported ML branch estimates. Record how many were adjusted.
  adjusted <- sum(initial$edge.length < 1e-8)
  initial$edge.length <- pmax(initial$edge.length, 1e-8)
  fit <- phangorn::pml(initial, data=data, k=gamma_categories, model=model)
  fit <- phangorn::optim.pml(fit, model=model, optNni=TRUE, optGamma=gamma_categories > 1L,
                              control=phangorn::pml.control(trace=0))
  if (!is.finite(fit$logLik) || any(!is.finite(fit$tree$edge.length)) || any(fit$tree$edge.length < 0)) {
    stop("ML fit did not produce finite, nonnegative branch estimates.")
  }
  replicates <- NULL
  if (bootstrap > 0L) {
    replicates <- phangorn::bootstrap.pml(fit, bs=bootstrap, optNni=TRUE, multicore=FALSE,
                                           control=phangorn::pml.control(trace=0))
    fit$tree$node.label <- as.character(round(100 * ape::prop.clades(fit$tree, replicates) / bootstrap, 1))
  }
  context <- fit$tree
  context$tip.label <- tip_data$TipLabel[match(context$tip.label, tip_data$InternalID)]
  names(alignment) <- tip_data$TipLabel
  asv_tree <- NULL
  if (length(queries) >= 2L) {
    asv_tree <- ape::keep.tip(fit$tree, internal[tip_data$Origin == "ASV"])
    asv_tree$tip.label <- tip_data$FeatureID[match(asv_tree$tip.label, tip_data$InternalID)]
    # Full-context clade support does not become ASV-only clade support on pruning.
    asv_tree$node.label <- NULL
  }
  if (!is.null(replicates)) {
    # bootstrap.pml may compress a multiPhylo's labels into a shared attribute.
    # Expand each tree before renaming, so that attribute cannot retain old IDs.
    replicates <- structure(lapply(replicates, function(tree) {
      tree$tip.label <- tip_data$TipLabel[match(tree$tip.label, tip_data$InternalID)]
      tree
    }), class="multiPhylo")
  }
  provenance <- list(method="joint_maximum_likelihood", model=model, gamma_categories=gamma_categories,
                      seed=seed, threads=threads, bootstrap=bootstrap, support_computed=bootstrap > 0,
                      min_overlap=min_overlap, starting_edges_adjusted=adjusted,
                      query_strands=query_strands, query_hashes=current_hashes,
                      reference_hashes=vapply(as.character(references), digest::digest, character(1), algo="xxhash64"),
                      matching_parameters=matches$parameters,
                      versions=c(DECIPHER=as.character(utils::packageVersion("DECIPHER")),
                                   phangorn=as.character(utils::packageVersion("phangorn"))))
  if (inherits(matches$reference, "silva_reference")) provenance$reference <- matches$reference
  structure(list(tree=context, asv_tree=asv_tree, alignment=alignment, tip_data=tip_data,
                   fit=fit, bootstrap_trees=replicates, provenance=provenance,
                   query_sequences=queries), class="asv_phylogeny")
}
