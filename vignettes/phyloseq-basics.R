## ---- eval=FALSE--------------------------------------------------------------
#  vignette("phyloseq_analysis")

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library("phyloseq")

## ---- eval=FALSE--------------------------------------------------------------
#  myOTU1 <- import_RDP_cluster("path/to/my/filename.clust")

## ---- eval=FALSE--------------------------------------------------------------
#  data(GlobalPatterns)
#  data(esophagus)
#  data(enterotype)
#  data(soilrep)

## -----------------------------------------------------------------------------
data(GlobalPatterns)
GlobalPatterns

## ---- eval=FALSE--------------------------------------------------------------
#  otu1 <- otu_table(raw_abundance_matrix, taxa_are_rows=FALSE)
#  sam1 <- sample_data(raw_sample_data.frame)
#  tax1 <- tax_table(raw_taxonomy_matrix)
#  tre1 <- read_tree(my_tree_file)

## ---- eval=FALSE--------------------------------------------------------------
#  ex1b <- phyloseq(my_otu_table, my_sample_data, my_taxonomyTable, my_tree)

## ---- eval=FALSE--------------------------------------------------------------
#  ex1c <- phyloseq(my_otu_table, my_sample_data)

## ----echo=FALSE---------------------------------------------------------------
topN <- 20

## -----------------------------------------------------------------------------
data(GlobalPatterns)
most_abundant_taxa <- sort(taxa_sums(GlobalPatterns), TRUE)[1:topN]
ex2 <- prune_taxa(names(most_abundant_taxa), GlobalPatterns)

## -----------------------------------------------------------------------------
topFamilies <- tax_table(ex2)[, "Family"]
as(topFamilies, "vector")

## ---- eval=FALSE--------------------------------------------------------------
#  testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
#  f1<- filterfun_sample(topk(2))
#  wh1 <- genefilter_sample(testOTU, f1, A=2)
#  wh2 <- c(T, T, T, F, F)
#  prune_taxa(wh1, testOTU)
#  prune_taxa(wh2, testOTU)

## -----------------------------------------------------------------------------
data(GlobalPatterns)
f1<- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(GlobalPatterns, f1, A=(1/2*nsamples(GlobalPatterns)))
sum(wh1)
ex2 <- prune_taxa(wh1, GlobalPatterns)

## -----------------------------------------------------------------------------
print(ex2)

## ---- eval=FALSE--------------------------------------------------------------
#  data(GlobalPatterns)
#  f1<- filterfun_sample(topf(0.9))
#  wh1 <- genefilter_sample(GlobalPatterns, f1, A=(1/3*nsamples(GlobalPatterns)))
#  sum(wh1)
#  prune_taxa(wh1, GlobalPatterns)

## -----------------------------------------------------------------------------
data("enterotype")
library("genefilter")
flist <- filterfun(kOverA(5, 2e-05))
ent.logi <- filter_taxa(enterotype, flist)
ent.trim <- filter_taxa(enterotype, flist, TRUE)
identical(ent.trim, prune_taxa(ent.logi, enterotype)) 
identical(sum(ent.logi), ntaxa(ent.trim))
filter_taxa(enterotype, flist, TRUE)

## -----------------------------------------------------------------------------
ex3 <- subset_samples(GlobalPatterns, SampleType%in%c("Freshwater", "Ocean", "Freshwater (creek)"))
ex3

## -----------------------------------------------------------------------------
subset(sample_data(GlobalPatterns), SampleType%in%c("Freshwater", "Ocean", "Freshwater (creek)"))

## -----------------------------------------------------------------------------
ex4 <- subset_taxa(GlobalPatterns, Phylum=="Firmicutes")
ex4

## -----------------------------------------------------------------------------
randomSpecies100 <- sample(taxa_names(GlobalPatterns), 100, replace=FALSE)
ex5 <- prune_taxa(randomSpecies100, GlobalPatterns)

## ---- eval=FALSE--------------------------------------------------------------
#  data(GlobalPatterns)
#  ex2 <- transform_sample_counts(GlobalPatterns, I)

## -----------------------------------------------------------------------------
ex4 <- transform_sample_counts(GlobalPatterns, threshrankfun(500))

## ---- eval=FALSE--------------------------------------------------------------
#  ex6 <- tax_glom(GlobalPatterns, taxrank = "Genus")

## ---- eval=FALSE--------------------------------------------------------------
#  ex7 <- tip_glom(GlobalPatterns, speciationMinLength = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("doParallel")
#  install.packages("doMC")
#  install.packages("doSNOW")
#  install.packages("doMPI")

