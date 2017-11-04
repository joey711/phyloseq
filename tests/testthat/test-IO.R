################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!

################################################################################
# import_mothur tests
mothlist  <- system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
mothgroup <- system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
mothtree  <- system.file("extdata", "esophagus.tree.gz", package="phyloseq")
cutoff    <- "0.10"
esophman  <- import_mothur(mothlist, mothgroup, mothtree, cutoff)	
# mothur "Shared" file, create with mothur from these example data files
mothshared = system.file("extdata", "esophagus.fn.shared.gz", package="phyloseq")
constaxonomy = system.file("extdata", "mothur_example.cons.taxonomy.gz", package="phyloseq")

test_that("import_mothur: import of esophagus dataset from mothur files in extdata/ produces a phyloseq object", {
	expect_that(esophman, is_a("phyloseq"))
})

test_that("import_mothur: The two phyloseq objects, example and just-imported, are identical", {
	data("esophagus")
	expect_that(esophagus, is_equivalent_to(esophman))
})

test_that("import_mothur: Test mothur file import on the (esophagus data).", {
	smlc <- show_mothur_cutoffs(mothlist)
	expect_that(smlc, is_equivalent_to(c("unique", "0.00", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.10")))	
})

test_that("import_mothur: abundances can be manipulated mathematically", {
	x1 <- as(otu_table(esophman), "matrix")
	expect_that(2*x1-x1, is_equivalent_to(x1) )
})

test_that("import_mothur: empty stuff is NULL", {
	expect_that(tax_table(esophman, FALSE), is_a("NULL"))
	expect_that(sample_data(esophman, FALSE), is_a("NULL"))
})

test_that("import_mothur: Expected classes of non-empty components", {
	expect_that(otu_table(esophman), is_a("otu_table"))
	expect_that(phy_tree(esophman), is_a("phylo"))
})

test_that("import_mothur: imported files become S4 object", {
	expect_that(isS4(esophman), is_true())
})

test_that("import_mothur: show method output tests", {
  expect_output(print(esophman), "phyloseq-class experiment-level object")
})

test_that("Test newer supported mothur output files, constaxonomy and shared files", {
  # Shared file
  sharedOTU = import_mothur(mothur_shared_file=mothshared, cutoff=cutoff)
  expect_is(sharedOTU, "otu_table")
  expect_equivalent(as(sharedOTU[1:5, ], "matrix")[, "B"], c(50, 37, 14, 52, 2))
  expect_equivalent(as(sharedOTU[1, ], "matrix")[1, ], c(50, 19, 5))
  # Constaxonomy file
  tax = import_mothur(mothur_constaxonomy_file=constaxonomy)
  expect_is(tax, "taxonomyTable")
})

################################################################################
# import_RDP tests
test_that("the import_RDP_otu function can properly read gzipped-example", {
  # Setup data
	otufile <- system.file("extdata",
	                       "rformat_dist_0.03.txt.gz",
	                       package="phyloseq")
	ex_otu  <- import_RDP_otu(otufile)	
	# test expectations
	expect_output(print(head(t(ex_otu))), "OTU Table:")
	expect_is(ex_otu, "otu_table")
	expect_equal(ntaxa(ex_otu), 5276)
	expect_equal(nsamples(ex_otu), 14)
	expect_is(sample_sums(ex_otu), "numeric")
})


################################################################################
# import_qiime tests
################################################################################
otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
mapfile <- system.file("extdata", "master_map.txt", package="phyloseq")
trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
rs_file <- system.file("extdata", "qiime500-refseq.fasta", package="phyloseq")

t0 <- import_qiime(otufile, mapfile, trefile, rs_file, verbose=FALSE)
test_that("Class of import result is phyloseq-class", {
	expect_that(t0, is_a("phyloseq"))
})

test_that("Classes of components are as expected", {
	expect_that(otu_table(t0), is_a("otu_table"))
	expect_that(tax_table(t0), is_a("taxonomyTable"))
	expect_that(sample_data(t0), is_a("sample_data"))
	expect_that(phy_tree(t0), is_a("phylo"))		
	expect_that(refseq(t0), is_a("DNAStringSet"))
})

test_that("Features of the abundance data are consistent, match known values", {
	expect_that(sum(taxa_sums(t0)), equals(1269671L))
	expect_that(sum(taxa_sums(t0)==0), equals(5L))
	expect_that(sum(taxa_sums(t0)>=100), equals(183L))
	expect_that(sum(taxa_sums(t0)), equals(sum(sample_sums(t0))))
	expect_that(sum(sample_sums(t0) > 10000L), equals(20L))
	expect_that(nsamples(t0), equals(26L))
	expect_that(ntaxa(t0), equals(500L))
	expect_that(length(rank_names(t0)), equals(7L))
})

test_that("Features of the taxonomy table match expected values", {
	expect_that(length(rank_names(t0)), equals(7L))
	expect_equal(rank_names(t0), 
               c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
	tax53 = as(tax_table(t0), "matrix")[53, ]
	expect_that(tax53, is_equivalent_to(c("Bacteria", "Proteobacteria", "Deltaproteobacteria",
		"Desulfovibrionales", "Desulfomicrobiaceae", 
    "Desulfomicrobium", "Desulfomicrobiumorale")))	
})
################################################################################
# parse function tests - note, these are also used by import_biom

test_that("Taxonomy vector parsing functions behave as expected", {
		
	chvec1 = c("Bacteria", "Proteobacteria", "Gammaproteobacteria",
		"Enterobacteriales", "Enterobacteriaceae", "Escherichia")
	
	chvec2 = c("k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
		"o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__")
	
	chvec3 = c("Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli",
		"o__Bacillales", "f__Staphylococcaceae")
	
	# Example where only some entries have greengenes prefix.
	chvec4 = c("Root", "k__Bacteria", "Firmicutes", "c__Bacilli",
		"o__Bacillales", "Staphylococcaceae", "z__mistake")

	# Even more terrible example, where leading or trailing space characters included
	# (the exact weirdnes of chvec4, compounded by leading and/or trailing space characters)
	chvec5 = c("  Root \n ", " k__Bacteria", "  Firmicutes", " c__Bacilli   ",
		"o__Bacillales  ", "Staphylococcaceae ", "\t z__mistake \t\n")		

	# This should give a warning because there were no greengenes prefixes
	expect_warning(t1 <- parse_taxonomy_greengenes(chvec1))
	# And output from previous call, t1, should be identical to default
	expect_that(parse_taxonomy_default(chvec1), is_equivalent_to(t1))
	
	# All the greengenes entries get trimmed by parse_taxonomy_greengenes
	expect_that(all(sapply(chvec2, nchar) > sapply(parse_taxonomy_greengenes(chvec2), nchar)), is_true())
	# None of the greengenes entries are trimmed by parse_taxonomy_default
	expect_that(any(sapply(chvec2, nchar) > sapply(parse_taxonomy_default(chvec2), nchar)), is_false())
	
	# Check that the "Root" element is not removed by parse_taxonomy_greengenes and parse_taxonomy_default.
	expect_that("Root" %in% chvec3, is_true())
	expect_that("Root" %in% parse_taxonomy_default(chvec3), is_true())
	expect_that(length(parse_taxonomy_default(chvec3)) == length(chvec3), is_true())
	
	# Check that non-greengenes prefixes, and those w/o prefixes, are given dummy rank(s)
	chvec4ranks = names(parse_taxonomy_greengenes(chvec4))
	expect_that(grep("Rank", chvec4ranks, fixed=TRUE), is_equivalent_to(c(1, 3, 6, 7)))
	# Check that everything given dummy rank in default parse.
	chvec4ranks = names(parse_taxonomy_default(chvec4))
	expect_that(grep("Rank", chvec4ranks, fixed=TRUE), is_equivalent_to(1:7))
	
	# chvec4 and chvec5 result in identical vectors.
	expect_that(parse_taxonomy_default(chvec4), is_equivalent_to(parse_taxonomy_default(chvec5)))
	expect_that(parse_taxonomy_greengenes(chvec4), is_equivalent_to(parse_taxonomy_greengenes(chvec5)))	
	
	# The names of chvec5, greengenes parsed, should be...
	correct5names = c("Rank1", "Kingdom", "Rank3", "Class", "Order", "Rank6", "Rank7")
	expect_that(names(parse_taxonomy_greengenes(chvec5)), is_equivalent_to(correct5names))
})

################################################################################
# import_biom tests

rich_dense_biom  <- system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
rich_sparse_biom <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
min_dense_biom   <- system.file("extdata", "min_dense_otu_table.biom",   package="phyloseq")
min_sparse_biom  <- system.file("extdata", "min_sparse_otu_table.biom",  package="phyloseq")
# the tree and refseq file paths that are suitable for all biom format style examples
treefilename = system.file("extdata", "biom-tree.phy",  package="phyloseq")
refseqfilename = system.file("extdata", "biom-refseq.fasta",  package="phyloseq")

test_that("Importing biom files yield phyloseq objects", {
	library("biomformat")
	rdbiom = read_biom(rich_sparse_biom)
	rsbiom = read_biom(rich_sparse_biom)

	rich_dense = import_biom(rdbiom)
	rich_sparse = import_biom(rsbiom)

	expect_that(rich_dense,  is_a("phyloseq"))
	expect_that(rich_sparse, is_a("phyloseq"))
	
	expect_that(ntaxa(rich_dense), equals(5L))
	expect_that(ntaxa(rich_sparse), equals(5L))

	# # Component classes
	# sample_data
	expect_that(access(rich_dense,  "sam_data"), is_a("sample_data"))
	expect_that(access(rich_sparse, "sam_data"), is_a("sample_data"))

	# taxonomyTable
	expect_that(access(rich_dense,  "tax_table"), is_a("taxonomyTable"))
	expect_that(access(rich_sparse, "tax_table"), is_a("taxonomyTable"))

	# otu_table
	expect_that(access(rich_dense,  "otu_table"), is_a("otu_table"))
	expect_that(access(rich_sparse, "otu_table"), is_a("otu_table"))
})

test_that("The different types of biom files yield phyloseq objects",{
	rich_dense = import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
	rich_sparse = import_biom(rich_sparse_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
	min_dense = import_biom(min_dense_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
	min_sparse = import_biom(min_sparse_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
	
	expect_that(rich_dense,  is_a("phyloseq"))
	expect_that(rich_sparse, is_a("phyloseq"))
	expect_that(min_dense,   is_a("phyloseq"))
	expect_that(min_sparse,  is_a("phyloseq"))

	expect_that(ntaxa(rich_dense), equals(5L))
	expect_that(ntaxa(rich_sparse), equals(5L))
	expect_that(ntaxa(min_dense), equals(5L))
	expect_that(ntaxa(min_sparse), equals(5L))			

	# # Component classes
	# sample_data
	expect_that(access(rich_dense,  "sam_data"), is_a("sample_data"))
	expect_that(access(rich_sparse, "sam_data"), is_a("sample_data"))
	expect_that(access(min_dense,   "sam_data"), is_a("NULL"))
	expect_that(access(min_sparse,  "sam_data"), is_a("NULL"))

	# taxonomyTable
	expect_that(access(rich_dense,  "tax_table"), is_a("taxonomyTable"))
	expect_that(access(rich_sparse, "tax_table"), is_a("taxonomyTable"))
	expect_that(access(min_dense,   "tax_table"), is_a("NULL"))
	expect_that(access(min_sparse,  "tax_table"), is_a("NULL"))		
	
	# phylo tree
	expect_that(access(rich_dense,  "phy_tree"), is_a("phylo"))
	expect_that(access(rich_sparse, "phy_tree"), is_a("phylo"))
	expect_that(access(min_dense,   "phy_tree"), is_a("phylo"))
	expect_that(access(min_sparse,  "phy_tree"), is_a("phylo"))

	# reference sequences
	expect_that(inherits(access(rich_dense,  "refseq"), "XStringSet"), is_true())
	expect_that(inherits(access(rich_sparse,  "refseq"), "XStringSet"), is_true())
	expect_that(inherits(access(min_dense,  "refseq"), "XStringSet"), is_true())
	expect_that(inherits(access(min_sparse,  "refseq"), "XStringSet"), is_true())
	expect_that(access(rich_dense,  "refseq"), is_a("DNAStringSet"))
	expect_that(access(rich_sparse, "refseq"), is_a("DNAStringSet"))
	expect_that(access(min_dense,   "refseq"), is_a("DNAStringSet"))
	expect_that(access(min_sparse,  "refseq"), is_a("DNAStringSet"))	
		
	# otu_table		
	expect_that(access(rich_dense,  "otu_table"), is_a("otu_table"))
	expect_that(access(rich_sparse, "otu_table"), is_a("otu_table"))
	expect_that(access(min_dense,   "otu_table"), is_a("otu_table"))
	expect_that(access(min_sparse,  "otu_table"), is_a("otu_table"))
	
	# Compare values in the otu_table. For some reason the otu_tables are not identical
	# one position is plus-two, another is minus-two
	combrich <- c(access(rich_dense, "otu_table"), access(rich_sparse, "otu_table"))
	expect_that(sum(diff(combrich, length(access(rich_dense, "otu_table")))), is_equivalent_to(0))
	expect_that(max(diff(combrich, length(access(rich_dense, "otu_table")))), is_equivalent_to(2))
	expect_that(min(diff(combrich, length(access(rich_dense, "otu_table")))), is_equivalent_to(-2))
	combmin <- c(access(min_dense, "otu_table"), access(min_sparse, "otu_table"))
	expect_that(sum(diff(combmin, length(access(min_dense, "otu_table")))), is_equivalent_to(0))
	expect_that(max(diff(combmin, length(access(min_dense, "otu_table")))), is_equivalent_to(2))
	expect_that(min(diff(combmin, length(access(min_dense, "otu_table")))), is_equivalent_to(-2))

	expect_that(access(min_dense, "otu_table"),  is_equivalent_to(access(rich_dense, "otu_table")))
	expect_that(access(min_sparse, "otu_table"), is_equivalent_to(access(rich_sparse, "otu_table")))

	# Compare values in the sample_data
	expect_that(access(rich_dense, "sam_data"), is_equivalent_to(access(rich_sparse, "sam_data")))
	
	# Compare values in the taxonomyTable
	expect_that(access(rich_dense, "tax_table"), is_equivalent_to(access(rich_sparse, "tax_table")))

})

test_that("the import_biom and import(\"biom\", ) syntax give same result", {
	x1 <- import_biom(rich_dense_biom, parseFunction=parse_taxonomy_greengenes)
	x2 <- import("biom", BIOMfilename=rich_dense_biom, parseFunction=parse_taxonomy_greengenes)	
	expect_that(x1, is_equivalent_to(x2))
})
################################################################################
# read_tree tests
test_that("The read_tree function works as expected:", {
	GPNewick <- read_tree(system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq"))
	expect_that(GPNewick, is_a("phylo"))
	expect_that(ntaxa(GPNewick), equals(length(GPNewick$tip.label)))
	expect_that(ntaxa(GPNewick), equals(500))
	expect_that(GPNewick$Nnode, equals(499))
	expect_that(taxa_names(GPNewick), is_equivalent_to(GPNewick$tip.label))	
	# Now read a nexus tree... 
	# Some error-handling expectations
	expect_that(read_tree("alskflsakjsfskfhas.akshfaksj"), gives_warning()) # file not exist
	not_tree <- system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
	expect_that(read_tree(not_tree), is_a("NULL")) # file not a tree, gives NULL
	expect_that(read_tree(not_tree, TRUE), throws_error()) # file not a tree, check turned on/TRUE
})
# read_tree_greengenes
test_that("The specialized read_tree_greengenes function works:", {
  # The included, gzipped version of the the 13_5 73% similarity greengenes tree.
  # It causes ape::read.tree to fail with an error, but read_tree_greengenes should be fine.
  treefile = system.file("extdata", "gg13-5-73.tree.gz", package="phyloseq")
  x = read_tree_greengenes(treefile)
  expect_is(x, "phylo")
  # Happen to know that all OTU names should be numbers.
  expect_match(x$tip.label, "[[:digit:]]+")
  # All tip/OTU names should be unique
  expect_false(any(duplicated(taxa_names(x))))
  # The more general read_tree function reads, but they are not equal
  y = read_tree(treefile)
  expect_false(all.equal(x, y))
})
################################################################################
# microbio_me_qiime tests
# This tests different features and expected behavior for
# the functioning of an interface function to the 
# microbio.me/qiime data repository.
#
zipfile = "study_816_split_library_seqs_and_mapping.zip"
zipfile = system.file("extdata", zipfile, package="phyloseq")
tarfile = "study_816_split_library_seqs_and_mapping.tar.gz"
tarfile = system.file("extdata", tarfile, package="phyloseq")
tarps = suppressWarnings(microbio_me_qiime(tarfile))
zipps = suppressWarnings(microbio_me_qiime(zipfile))
# This function is intended to interface with an external server,
# as described in the documentation.
# However, I don't want successful testing of this package to 
# rely on the presence or form of particular files on an 
# external server, so these tests will be done exclusively on 
# compressed file(s) representing what is exposed by the data server
# It is up to the user to provide valid a URL in practice,
# and the function attempts to provide informative status
# and error messages if things go awry.
test_that("The microbio_me_qiime imports as expected: .tar.gz", {
  expect_is(tarps, "phyloseq")
  expect_is(sample_data(tarps, errorIfNULL=FALSE), "sample_data")
  expect_is(otu_table(tarps, errorIfNULL=FALSE), "otu_table")
  expect_identical(nrow(otu_table(tarps)), 50L)
  expect_identical(nrow(sample_data(tarps)), 15L)
})
test_that("The microbio_me_qiime imports as expected: .zip", {  
  expect_is(zipps, "phyloseq")
  expect_is(sample_data(zipps, errorIfNULL=FALSE), "sample_data")
  expect_is(otu_table(zipps, errorIfNULL=FALSE), "otu_table")
  expect_identical(nrow(otu_table(zipps)), 50L)
  expect_identical(nrow(sample_data(zipps)), 15L)
})
test_that("Results of .tar.gz and .zip should be identical", {  
  expect_identical(tarps, zipps)
  expect_identical(sample_data(tarps), sample_data(zipps))
  expect_identical(otu_table(tarps), otu_table(zipps))
})
################################################################################
# import_usearch_uc
################################################################################
usearchfile = system.file("extdata", "usearch.uc", package="phyloseq")
OTU1 = import_usearch_uc(usearchfile)
test_that("import_usearch_uc: Properly omit entries from failed search", {  
  ucLines = readLines(usearchfile)
  expect_identical( sum(OTU1), (length(ucLines) - length(grep("*", ucLines, fixed=TRUE))) )
  expect_identical( nrow(OTU1), 37L)
  expect_identical( nrow(OTU1), nsamples(OTU1))
  expect_identical( ncol(OTU1), ntaxa(OTU1))
  expect_identical( ncol(OTU1), 33L)
  expect_equivalent(colSums(OTU1)[1:6], c(6, 1, 2, 1, 1, 1))
})
################################################################################