library("phyloseq")
library("testthat")

# import biom test
silva_biom <- system.file("extdata", "SILVA_OTU_table.biom", package="phyloseq")
# Create phyloseq object using silva parseing function
silva_phyloseq <- import_biom(BIOMfilename = silva_biom, parseFunction = parse_taxonomy_silva_128)
