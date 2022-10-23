## ----load-phyloseq, message=FALSE, warning=FALSE------------------------------
library("phyloseq"); packageVersion("phyloseq")

## ----filepath-----------------------------------------------------------------
filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
kostic = microbio_me_qiime(filepath)

## ----example-path-local, eval=FALSE-------------------------------------------
#  filepath = "~/Downloads/study_1457_split_library_seqs_and_mapping.zip"
#  kostic = microbio_me_qiime(filepath)

## ----example-path-remote, eval=FALSE------------------------------------------
#  kostic = microbio_me_qiime(1457)

## ----show-variables-----------------------------------------------------------
kostic
head(sample_data(kostic)$DIAGNOSIS, 10)

## ----deseq2, message=FALSE, warning=FALSE-------------------------------------
library("DESeq2"); packageVersion("DESeq2")

## ----rm-bad-samples-----------------------------------------------------------
kostic <- subset_samples(kostic, DIAGNOSIS != "None")
kostic <- prune_samples(sample_sums(kostic) > 500, kostic)
kostic

## ----run-deseq2---------------------------------------------------------------
diagdds = phyloseq_to_deseq2(kostic, ~ DIAGNOSIS)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

## ----grab-results-process-table-----------------------------------------------
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)

## ----table-prelim-------------------------------------------------------------
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

## ----make-markdown-table, echo=FALSE, results='asis'--------------------------
# Make a markdown table
posigtab = data.frame(OTU=rownames(posigtab), posigtab)
cat(paste(colnames(posigtab), collapse=" | "), fill=TRUE)
cat(paste(rep("---", times=ncol(posigtab)), collapse=" | "), fill=TRUE)
dummy = apply(posigtab, 1, function(x){
  cat(paste(x, collapse=" | "), fill=TRUE)
})

## ----bar-plot-----------------------------------------------------------------
library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

