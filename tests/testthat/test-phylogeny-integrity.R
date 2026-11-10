# Independent reference using root paths of individual tips, rather than the
# internal-node aggregation used by fastUniFrac.
reference_unifrac <- function(tree, counts, weighted=FALSE, normalized=TRUE) {
  incidence <- matrix(0, nrow(tree$edge), nrow(counts),
                       dimnames=list(NULL, rownames(counts)))
  for (tip in seq_along(tree$tip.label)) {
    current <- tip
    repeat {
      edge <- which(tree$edge[, 2] == current)
      if (!length(edge)) break
      incidence[edge, tree$tip.label[tip]] <- 1
      current <- tree$edge[edge, 1]
    }
  }
  mass <- incidence %*% sweep(counts, 2, colSums(counts), "/")
  answer <- matrix(0, ncol(counts), ncol(counts), dimnames=list(colnames(counts), colnames(counts)))
  for (a in seq_len(ncol(counts) - 1L)) for (b in seq.int(a + 1L, ncol(counts))) {
    numerator <- if (weighted) abs(mass[, a] - mass[, b]) else xor(mass[, a] > 0, mass[, b] > 0)
    denominator <- if (weighted) mass[, a] + mass[, b] else mass[, a] > 0 | mass[, b] > 0
    value <- sum(tree$edge.length * numerator)
    if (!weighted || normalized) value <- value / sum(tree$edge.length * denominator)
    answer[a, b] <- answer[b, a] <- value
  }
  stats::as.dist(answer)
}

test_that("UniFrac handles polytomies and unary nodes without changing topology", {
  for (newick in c("((a:1,b:1,c:1):1,d:2);", "(((a:1):0.5,b:1):1,c:2);",
                   "((a:0.7,b:2):0.5,(c:1,d:0.2):1.3);")) {
    tree <- ape::read.tree(text=newick)
    values <- matrix(seq_len(length(tree$tip.label) * 3) %% 5, length(tree$tip.label),
                       dimnames=list(tree$tip.label, c("01", "u", "v")))
    values <- values[rev(rownames(values)), , drop=FALSE]
    for (rows in c(TRUE, FALSE)) {
      input <- phyloseq(otu_table(if (rows) values else t(values), rows), tree)
      original <- input
      for (weighted in c(TRUE, FALSE)) for (normalized in c(TRUE, FALSE)) {
        expected <- reference_unifrac(tree, values, weighted, normalized)
        expect_warning(actual <- UniFrac(input, weighted=weighted, normalized=normalized), NA)
        expect_equal(as.matrix(actual), as.matrix(expected), tolerance=1e-12)
      }
      expect_identical(input, original)
      expect_identical(phy_tree(input), tree)
    }
  }
})

test_that("polytomy regression agrees with exact answers and zero-edge resolutions", {
  tree <- ape::read.tree(text="((a:1,b:1,c:1):1,d:2);")
  values <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0), 4,
                     dimnames=list(c("a", "b", "c", "d"), c("u", "v")))
  input <- phyloseq(otu_table(values, TRUE), tree)
  resolved <- phyloseq(otu_table(values, TRUE), ape::multi2di(tree, random=FALSE))
  expect_equal(as.numeric(UniFrac(input)), 2/3)
  expect_equal(as.numeric(UniFrac(input, weighted=TRUE)), 0.5)
  expect_equal(as.numeric(UniFrac(input, weighted=TRUE, normalized=FALSE)), 2)
  for (weighted in c(TRUE, FALSE)) for (normalized in c(TRUE, FALSE)) {
    expect_equal(UniFrac(input, weighted=weighted, normalized=normalized),
                   UniFrac(resolved, weighted=weighted, normalized=normalized))
  }
})

test_that("merging incomplete refseq rejects silent taxon loss or explicitly drops refseq", {
  for (rows in c(TRUE, FALSE)) {
    a <- matrix(10, 1, 1, dimnames=list("asv1", "a"))
    b <- matrix(c(10, 4), 2, 1, dimnames=list(c("asv1", "asv2"), "b"))
    meta <- function(id) sample_data(data.frame(group="case", batch=1, row.names=id))
    x <- phyloseq(otu_table(if (rows) a else t(a), rows), meta("a"),
                   Biostrings::DNAStringSet(c(asv1="ACGT")))
    y <- phyloseq(otu_table(if (rows) b else t(b), rows), meta("b"))
    original <- x
    expect_error(merge_phyloseq(x, y), "missing for 1 taxa.*asv2.*refseq_policy")
    expect_error(merge_phyloseq(y, x), "missing for 1 taxa")
    expect_warning(merged <- merge_phyloseq(x, y, refseq_policy="drop"), "Dropping.*refseq")
    expect_true(methods::validObject(merged))
    expect_equal(sum(sample_sums(merged)), 24)
    expect_setequal(taxa_names(merged), c("asv1", "asv2"))
    expect_null(refseq(merged, errorIfNULL=FALSE))
    expect_identical(taxa_are_rows(merged), rows)
    expect_identical(x, original)
    # Complete sequence unions remain supported with the normal default.
    sequences <- Biostrings::DNAStringSet(c(asv1="ACGT", asv2="AAAA"))
    full <- phyloseq(otu_table(y), sample_data(y), sequences)
    expect_warning(complete <- merge_phyloseq(x, full), NA)
    expect_equal(sum(sample_sums(complete)), 24)
    expect_setequal(names(refseq(complete)), c("asv1", "asv2"))
  }
})
