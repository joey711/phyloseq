# Opt-in integration check with real SILVA files, not a biological accuracy test.
# Rscript --vanilla .github/diagnostics/validate-silva-live.R CACHE OUTPUT [--download]
# The explicit flag permits large downloads; otherwise the cache must exist.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 3L ||
    (length(args) == 3L && args[[3L]] != "--download")) {
  stop("Usage: validate-silva-live.R CACHE OUTPUT [--download]")
}
cache <- path.expand(args[[1L]])
output <- path.expand(args[[2L]])
dir.create(output,recursive=TRUE,showWarnings=FALSE)
pkgload::load_all(".",quiet=TRUE)
timings <- list()
timings$reference <- system.time(raw <- silva_reference("144",cache_dir=cache,
                                                       download=length(args)==3L))
saveRDS(raw,file.path(output,"source-reference.rds"))
cat("Publisher checksum validation passed.\n")
flush.console()
qc_directory <- file.path(output,"qc")
manifest <- file.path(qc_directory,"reference.rds")
if (file.exists(manifest)) {
  filtered <- readRDS(manifest)
  stopifnot(identical(filtered$source_reference$checksums,raw$checksums))
} else {
  timings$qc <- system.time(filtered <- filter_silva_reference(raw,qc_directory))
}
print(filtered$summary)
flush.console()
# Test the complete on-disk database search with controlled, known-origin
# fragments from the first unambiguous passing genes. This is a plumbing check:
# fragments are reference-derived, not independently measured biological ASVs.
references <- Biostrings::readDNAStringSet(filtered$fasta,nrec=100)
names(references) <- sub("[[:space:]].*$","",names(references))
references <- references[!grepl("[^ACGT]",as.character(references))]
stopifnot(length(references)>=3L)
parents <- references[rep(1L,2L)]
# Shared coordinates on the same parent guarantee known biological overlap.
# Equal raw offsets on different partial genes would not guarantee homology.
queries <- Biostrings::subseq(parents,start=c(300L,325L),width=250L)
names(queries) <- c("derived_01","derived_02")
timings$matching <- system.time(hits <- match_asvs_to_reference(queries,filtered,max_hits=100))
stopifnot(all(hits$summary$Status=="matched"),
          all(names(parents) %in% hits$matches$ReferenceID | hits$summary$Truncated))
write.table(hits$summary,file.path(output,"query-summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(hits$matches,file.path(output,"candidates.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("Full-database search passed.\n")
flush.console()
# Keep inference targeted. Use three or more actually retrieved references and
# declare the context explicitly when candidate caps were reached.
chosen <- unique(unlist(lapply(names(queries),function(id) {
  head(hits$matches$ReferenceID[hits$matches$QueryID==id],10L)
}),use.names=FALSE))
if (length(chosen)>=3L && all(hits$summary$Strand %in% c("+","-")) &&
    all(names(queries) %in% hits$matches$QueryID[hits$matches$ReferenceID %in% chosen])) {
  timings$inference <- system.time(result <- infer_asv_tree(hits,reference_ids=chosen,
                                                            bootstrap=0L,seed=42L))
  stopifnot(setequal(result$asv_tree$tip.label,names(queries)),
            identical(as.character(result$query_sequences),as.character(queries)))
  saveRDS(result,file.path(output,"phylogeny.rds"))
  Biostrings::writeXStringSet(result$alignment,file.path(output,"alignment.fasta"))
  cat("Targeted joint ML inference passed (support not assessed).\n")
} else {
  stop("Search succeeded, but this reference-derived fixture needs an explicit context/orientation decision before inference.")
}
report <- list(release=raw$version,source_checksums=raw$checksums,
                 qc=filtered$summary,query_summary=hits$summary,
                 reference_ids=chosen,timings=timings,
                 inference=result$provenance,session=capture.output(sessionInfo()),
                 scope="Reference-derived fragments: integration plumbing only; not biological accuracy.")
saveRDS(report,file.path(output,"validation.rds"))
print(timings)
writeLines(report$session,file.path(output,"session-info.txt"))
cat("Live integration validation saved to:",normalizePath(output),"\n")
