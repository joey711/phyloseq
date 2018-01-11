library("phyloseq")
library("testthat")

# import biom test
silva_biom <- system.file("extdata", "SILVA_OTU_table.biom", package="phyloseq")
# Create phyloseq object using silva parseing function
silva_phyloseq <- import_biom(BIOMfilename = silva_biom, parseFunction = parse_taxonomy_silva_128)

test_that("parse_taxonomy_silva_128: silva parser produces a phyloseq object.", {
  expect_that(silva_phyloseq, is_a("phyloseq"))
})

test_that("parse_taxonomy_silva_128: imported files become S4 object", {
  expect_that(isS4(silva_phyloseq), is_true())
})

test_that("Classes of components are as expected", {
  expect_that(otu_table(silva_phyloseq), is_a("otu_table"))
  expect_that(tax_table(silva_phyloseq), is_a("taxonomyTable"))
})

test_that("Features of the taxonomy table match expected values", {
  expect_that(length(rank_names(silva_phyloseq)), equals(7L))
  expect_equal(rank_names(silva_phyloseq), 
               c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  tax54 = as(tax_table(silva_phyloseq), "matrix")[54, ]
  expect_that(tax54, is_equivalent_to(c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae","Prevotella 9", "uncultured bacterium")))	
})