silva_fixture <- function(directory, version="144") {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  repeat_dna <- function(n) substr(paste(rep("ACGT", ceiling(n/4)), collapse=""), 1, n)
  ids <- c("001.1.1200", "arch.1.900", "short.1.1199", "euk.1.1800", "amb.1.1200",
            "seq.1.1200", "ali.1.1200", "hom.1.1200", "vec.1.1200", "missing.1.1200",
            "boundary.1.1200", "001.1301.2500")
  dna <- c(repeat_dna(1200), repeat_dna(900), repeat_dna(1199), repeat_dna(1800),
             paste0(repeat_dna(1170), paste(rep("N", 30), collapse="")), rep(repeat_dna(1200), 5),
             paste0(repeat_dna(1176), paste(rep("N", 24), collapse="")), repeat_dna(1200))
  domain <- c("Bacteria", "Archaea", "Bacteria", "Eukaryota", rep("Bacteria", 8))
  names(dna) <- paste(ids, paste0(domain, ";Test"))
  fasta <- file.path(directory, "reference.fasta")
  Biostrings::writeXStringSet(Biostrings::BStringSet(dna), fasta)
  parts <- strsplit(ids, ".", fixed=TRUE)
  quality <- data.frame("Primary Accession"=vapply(parts, `[`, character(1), 1),
                         Start=vapply(parts, `[`, character(1), 2), Stop=vapply(parts, `[`, character(1), 3),
                         "Sequence Quality"=100, "Alignment Quality"=100,
                         "% Ambiguities"=0, "% Homopolymers"=0, "% Vector Contamination"=0,
                         check.names=FALSE)
  quality[["Sequence Quality"]][6] <- 94
  quality[["Alignment Quality"]][7] <- 74
  quality[["% Homopolymers"]][8] <- 3
  quality[["% Vector Contamination"]][9] <- 3
  quality[["Sequence Quality"]][10] <- NA
  quality[["% Ambiguities"]][11] <- 2
  if (version == "138.2") quality[["Pintail Quality"]] <- c(90, rep(100, 11))
  file <- file.path(directory, "quality.tsv")
  utils::write.table(quality, file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  list(reference=silva_reference(version, fasta=fasta, quality=file), ids=ids,
         quality=quality, fasta=fasta, quality_file=file)
}

asv_fixture <- function() {
  values <- matrix(1:9, 3, dimnames=list(c("01", "asv2", "asv3"), c("01", "s2", "s3")))
  phyloseq(otu_table(values, TRUE),
             sample_data(data.frame(group=c("a", "b", "a"), row.names=colnames(values))),
             Biostrings::DNAStringSet(c("01"="ACGTAC",asv2="TTAAAC",asv3="ACGTAC")))
}

test_that("ASV selection retains duplicated DNA and aligns counts, sequences and samples", {
  input <- asv_fixture()
  original <- input
  by_sequence <- select_asvs(input, sequences="acguac")
  expect_identical(taxa_names(by_sequence), c("01", "asv3"))
  expect_identical(as(otu_table(by_sequence), "matrix"), as(otu_table(input), "matrix")[c(1,3), , drop=FALSE])
  expect_identical(sample_data(by_sequence), sample_data(input))
  expect_identical(names(refseq(by_sequence)), c("01", "asv3"))
  expect_identical(select_asvs(input, ids=c("asv3", "01")), by_sequence)
  expect_identical(input, original)
})

test_that("ASV selectors validate IDs and orientation without silent loss", {
  input <- asv_fixture()
  original <- input
  expect_error(select_asvs(input, sequences="GTAAAC"), "do not match")
  expect_identical(taxa_names(select_asvs(input, sequences="GTTTAA", both_strands=TRUE)), "asv2")
  expect_identical(as.character(asv_sequences(c("01"="acgu"))), c("01"="ACGT"))
  expect_identical(names(asv_sequences(input, ids=c("asv3", "01"))), c("asv3", "01"))
  expect_error(select_asvs(input), "exactly one")
  expect_error(select_asvs(input, ids="01", sequences="ACGTAC"), "exactly one")
  expect_error(select_asvs(input, ids="missing"), "absent")
  expect_error(select_asvs(input, ids=c("01","01")), "unique")
  expect_error(select_asvs(input, ids="01", prune=NA), "prune")
  expect_error(asv_sequences(c(a="ACGN")), "unambiguous")
  expect_error(asv_sequences(c(a="AC-G")), "unambiguous")
  expect_error(asv_sequences(c(a="")), "unambiguous")
  expect_error(asv_sequences(c("ACGT")), "IDs")
  expect_error(asv_sequences(phyloseq(otu_table(input),sample_data(input))), "missing")
  expect_identical(input, original)
})

test_that("TSE ASV selection retains assays and stale-transformation diagnostics", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- transform_composition(asv_fixture(), "clr")
  original <- input
  filtered <- select_asvs(input, sequences="ACGTAC")
  expect_s4_class(filtered, "TreeSummarizedExperiment")
  expect_identical(rownames(filtered), c("01", "asv3"))
  expect_identical(as.character(asv_sequences(filtered)), as.character(asv_sequences(input, ids=c("01", "asv3"))))
  expect_identical(SummarizedExperiment::colData(filtered), SummarizedExperiment::colData(input))
  expect_false(check_composition(filtered, assay="clr", scale="clr")$valid)
  expect_identical(input, original)
})

test_that("16S QC uses domain-specific length minima and explicit score boundaries", {
  fixture <- silva_fixture(tempfile())
  output <- tempfile()
  filtered <- filter_silva_reference(fixture$reference, output, chunk_size=2)
  expect_s3_class(filtered, "silva_filtered_reference")
  expect_equal(filtered$summary$total, 12)
  expect_equal(filtered$summary$retained, 4)
  report <- utils::read.delim(filtered$qc_report, stringsAsFactors=FALSE)
  expected <- fixture$ids[c(1,2,11,12)]
  expect_identical(report$ReferenceID[report$Keep], expected)
  expect_identical(sub(" .*$","",names(Biostrings::readDNAStringSet(filtered$fasta))), expected)
  expect_equal(report$AmbiguousFraction[11], 0.02)
  expect_match(report$Reason[3], "length")
  expect_match(report$Reason[4], "domain")
  expect_match(report$Reason[5], "ambiguity")
  expect_match(report$Reason[6], "sequence_quality")
  expect_match(report$Reason[7], "alignment_quality")
  expect_match(report$Reason[8], "homopolymer")
  expect_match(report$Reason[9], "vector")
  expect_match(report$Reason[10], "missing_quality")
  expect_identical(unname(tools::md5sum(fixture$fasta)), fixture$reference$checksums$fasta)
  expect_false(filtered$publisher_verified)
  expect_identical(filtered$source_reference, fixture$reference)
  expect_error(filter_silva_reference(fixture$reference, output), "already exist")
  archaea <- filter_silva_reference(fixture$reference, tempfile(), domains="Archaea")
  expect_equal(archaea$summary$retained, 1)
  none <- filter_silva_reference(fixture$reference, tempfile(), min_sequence_quality=100,
                                  min_length=c(Bacteria=2000,Archaea=2000))
  expect_equal(none$summary$retained, 0)
})

test_that("quality schema, missing data and release differences are enforced", {
  fixture <- silva_fixture(tempfile())
  expect_error(filter_silva_reference(fixture$reference,tempfile(),min_pintail=90), "discontinued")
  expect_error(filter_silva_reference(fixture$reference,tempfile(),missing_quality="error"), "missing")
  older <- silva_fixture(tempfile(),version="138.2")
  filtered <- filter_silva_reference(older$reference,tempfile(),min_pintail=95)
  report <- utils::read.delim(filtered$qc_report)
  expect_match(report$Reason[1], "pintail")
  expect_false(report$Keep[1])
  expect_error(filter_silva_reference(fixture$reference,tempfile(),domains="Eukaryota"), "domains")
  expect_error(filter_silva_reference(fixture$reference,tempfile(),max_ambiguous=NA), "max_ambiguous")
  expect_error(filter_silva_reference(fixture$reference,tempfile(),min_length=1200), "min_length")
  expect_error(filter_silva_reference(fixture$reference,tempfile(),chunk_size=0), "chunk_size")
  cat("\nchanged",file=fixture$fasta,append=TRUE)
  expect_error(filter_silva_reference(fixture$reference,tempfile()), "changed|checksum")
  fixture <- silva_fixture(tempfile())
  altered <- fixture$quality
  altered[["Sequence Quality"]][1] <- -1
  utils::write.table(altered,fixture$quality_file,sep="\t",quote=FALSE,row.names=FALSE)
  fresh <- silva_reference("144",fasta=fixture$fasta,quality=fixture$quality_file)
  expect_error(filter_silva_reference(fresh,tempfile()), "0-100")
  altered <- fixture$quality[,names(fixture$quality) != "Alignment Quality"]
  utils::write.table(altered,fixture$quality_file,sep="\t",quote=FALSE,row.names=FALSE)
  fresh <- silva_reference("144",fasta=fixture$fasta,quality=fixture$quality_file)
  expect_error(filter_silva_reference(fresh,tempfile()), "unavailable")
})

test_that("reference cache checks digests and never downloads implicitly", {
  fixture <- silva_fixture(tempfile())
  expect_error(silva_reference("144",cache_dir=tempfile()), "not cached")
  expect_error(silva_reference("132",cache_dir=tempfile()), "Official download")
  expect_error(silva_reference("../144"), "release string")
  expect_error(silva_reference("144",fasta=fixture$fasta), "both local")
  cache <- tempfile()
  directory <- file.path(cache,"SILVA_144_SSURef_NR99")
  dir.create(directory,recursive=TRUE)
  files <- c(fasta="SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz",quality="SILVA_144_SSURef_Nr99.quality.gz")
  # Offline cache fixture: checksum mechanics are exercised without claiming that
  # these synthetic contents are an authentic SILVA release.
  for (field in names(files)) {
    destination <- file.path(directory,files[[field]])
    file.copy(fixture$reference[[field]],destination)
    writeLines(paste(unname(tools::md5sum(destination)),basename(destination)),paste0(destination,".md5"))
  }
  cached <- silva_reference("144",cache_dir=cache)
  expect_identical(cached$version,"144")
  expect_true(cached$publisher_verified)
  expect_match(cached$urls[["fasta"]],"release_144/Exports/SILVA_144")
  writeLines("bad",cached$fasta)
  expect_error(silva_reference("144",cache_dir=cache), "checksum")
})

test_that("compressed exports stream into QC and matching without losing reference IDs", {
  fixture <- silva_fixture(tempfile())
  fasta <- paste0(fixture$fasta,".gz")
  quality <- paste0(fixture$quality_file,".gz")
  Biostrings::writeXStringSet(Biostrings::readBStringSet(fixture$fasta),fasta,compress=TRUE)
  connection <- gzfile(quality,"wt")
  writeLines(readLines(fixture$quality_file),connection)
  close(connection)
  raw <- silva_reference("144",fasta=fasta,quality=quality)
  filtered <- filter_silva_reference(raw,tempfile(),chunk_size=1)
  candidates <- match_asvs_to_reference(c(q="ACGTACGT"),filtered,chunk_size=1,max_hits=Inf)
  expect_setequal(candidates$matches$ReferenceID,fixture$ids[c(1,2,11,12)])
  expect_false(candidates$summary$Truncated)
  expect_identical(candidates$summary$Strand,"ambiguous")
})

test_that("cross-chunk duplicate IDs fail QC without publishing incomplete outputs", {
  fixture <- silva_fixture(tempfile())
  first <- Biostrings::readBStringSet(fixture$fasta,nrec=1)
  Biostrings::writeXStringSet(first,fixture$fasta,append=TRUE)
  raw <- silva_reference("144",fasta=fixture$fasta,quality=fixture$quality_file)
  output <- tempfile()
  expect_error(filter_silva_reference(raw,output,chunk_size=1),"Duplicate reference")
  expect_length(list.files(output,all.files=TRUE,no..=TRUE),0)
})

test_that("candidate matching preserves ambiguity, strand and explicit truncation", {
  refs <- c(ref1="TTACGTACTT",ref2="GGACGTACAA",reverse="GGGTAAACGG")
  queries <- c("01"="ACGTAC",rev="GTTTAC",missing="AAAAAA")
  result <- match_asvs_to_reference(queries,refs)
  expect_equal(result$summary$Candidates,c(2,1,0))
  expect_identical(result$summary$Strand,c("+","-",NA_character_))
  expect_identical(result$summary$Status,c("matched","matched","no_hit"))
  expect_true(all(result$matches$Identity==1))
  expect_true(all(result$matches$QueryCoverage==1))
  expect_identical(as.character(result$query_sequences),queries)
  capped <- match_asvs_to_reference(queries,refs,max_hits=1)
  expect_true(capped$summary$Truncated[1])
  expect_equal(sum(capped$matches$QueryID=="01"),1)
  approximate <- match_asvs_to_reference(c(q="ACGTAT"),refs,max_mismatches=1,both_strands=FALSE)
  expect_true(all(approximate$matches$Mismatches==1))
  expect_equal(approximate$matches$Identity,rep(5/6,nrow(approximate$matches)))
  palindromic <- match_asvs_to_reference(c(q="ACGT"),c(ref="TTACGTTT"))
  expect_identical(palindromic$summary$Strand,"ambiguous")
  expect_error(match_asvs_to_reference(queries,refs,max_mismatches=10),"max_mismatches")
  fixture <- silva_fixture(tempfile())
  expect_error(match_asvs_to_reference(queries,fixture$reference),"Filter")
})

phylogeny_fixture <- function() {
  set.seed(414)
  bases <- sample(c("A","C","G","T"),160,TRUE)
  references <- vapply(1:4,function(i) {
    b <- bases
    b[c(i,130+i)] <- c("A","C","G","T")[i]
    paste(b,collapse="")
  },character(1))
  names(references) <- c("01","r2","r3","r4")
  queries <- c(substr(references[1],30,100),substr(references[3],50,120))
  names(queries) <- c("01","asv2")
  match_asvs_to_reference(queries,references)
}

test_that("joint inference preserves observed DNA, IDs, RNG and valid branch estimates", {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("phangorn")
  matches <- phylogeny_fixture()
  original <- matches
  rng <- .Random.seed
  fit <- infer_asv_tree(matches,model="JC",gamma_categories=1,bootstrap=2)
  expect_s3_class(fit,"asv_phylogeny")
  expect_identical(.Random.seed,rng)
  expect_identical(matches,original)
  expect_identical(fit$query_sequences,matches$query_sequences)
  expect_setequal(fit$tree$tip.label,c(paste0("REF::",names(matches$reference)),paste0("ASV::",names(matches$query_sequences))))
  expect_setequal(fit$asv_tree$tip.label,names(matches$query_sequences))
  expect_false(ape::is.rooted(fit$tree))
  expect_true(all(is.finite(fit$tree$edge.length) & fit$tree$edge.length>=0))
  ungapped <- gsub("[-.]","",as.character(fit$alignment))
  expect_identical(unname(ungapped[fit$tip_data$Origin=="ASV"]),unname(as.character(matches$query_sequences)))
  expect_equal(length(fit$bootstrap_trees),2)
  expect_true(all(vapply(fit$bootstrap_trees,function(tree) setequal(tree$tip.label,fit$tree$tip.label),logical(1))))
  expect_null(fit$asv_tree$node.label)
  expect_true(fit$provenance$support_computed)
  expect_identical(fit$provenance$method,"joint_maximum_likelihood")
})

test_that("joint inference rejects missing, truncated or altered evidence", {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("phangorn")
  matches <- phylogeny_fixture()
  expect_error(infer_asv_tree(matches,reference_ids="01"),"three reference")
  expect_error(infer_asv_tree(matches,max_sequences=3),"Too many")
  expect_error(infer_asv_tree(matches,query_strands=c("01"="-",asv2="-")),"chosen query orientation")
  expect_error(infer_asv_tree(matches,min_overlap=200),"Insufficient shared")
  capped <- match_asvs_to_reference(matches$query_sequences,matches$reference,max_hits=2)
  expect_error(infer_asv_tree(capped),"truncated")
  altered <- matches
  altered$matches$Mismatches[1] <- 10L
  expect_error(infer_asv_tree(altered),"Candidate matches")
  altered <- matches
  altered$reference[1] <- paste(rep("A",160),collapse="")
  expect_error(infer_asv_tree(altered),"Reference sequences changed")
  matches$query_sequences[[1]] <- Biostrings::DNAString(paste(rep("A",71),collapse=""))
  expect_error(infer_asv_tree(matches),"changed")
})

test_that("default GTR fit restores the absence of an RNG state", {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("phangorn")
  matches <- phylogeny_fixture()
  saved_rng <- .Random.seed
  on.exit(assign(".Random.seed",saved_rng,envir=.GlobalEnv),add=TRUE)
  rm(".Random.seed",envir=.GlobalEnv)
  fit <- infer_asv_tree(matches)
  expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
  expect_identical(fit$provenance$model,"GTR")
  expect_false(fit$provenance$support_computed)
  expect_true(is.finite(fit$fit$logLik))
})
