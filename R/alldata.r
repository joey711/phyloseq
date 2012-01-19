################################################################################
#' (Data) Example data set 1 for the phyloseq package. (2011)
#'
#' This is a small preliminary dataset of human gut microbial communities from
#' subjects of different gender and diet. The exact meaning of diet and gender
#' has been obscured because the larger study to which it belongs has not yet
#' been published. For a more complete and already-published example dataset,
#' try the \code{\link{enterotype}} data.
#'
#' @name data-ex1
#' @aliases ex1
#' @docType data
#' @author Paul J. McMurdie II \email{mcmurdie@@stanford.edu}
#' @references \url{www.stanford.edu/~mcmurdie}
#' @keywords data
#' @examples
#' ## data(ex1)
################################################################################
# This is a dummy line. This source file is just for documenting the
# example data, ex1, for the phyloseq package.
NA
################################################################################
#' (Data) Small example dataset from a human esophageal community (2004)
#'
#' Includes just 3 samples, 1 each from 3 subjects. Although the research article mentions 4 subjects,
#' only 3 are included in this dataset.
#'
#' abstract from research article (quoted):
#' 
#' The esophagus, like other luminal organs of the digestive system, provides a potential environment for bacterial colonization, but little is known about the presence of a bacterial biota or its nature. By using broad-range 16S rDNA PCR, biopsies were examined from the normal esophagus of four human adults. The 900 PCR products cloned represented 833 unique sequences belonging to 41 genera, or 95 species-level operational taxonomic units (SLOTU); 59 SLOTU were homologous with culture-defined bacterial species, 34 with 16S rDNA clones, and two were not homologous with any known bacterial 16S rDNA. Members of six phyla, Firmicutes, Bacteroides, Actinobacteria, Proteobacteria, Fusobacteria, and TM7, were represented. A large majority of clones belong to 13 of the 41 genera (783/900, 87\%), or 14 SLOTU (574/900, 64\%) that were shared by all four persons. Streptococcus (39\%), Prevotella (17\%), and Veilonella (14\%) were most prevalent. The present study identified 56-79\% of SLOTU in this bacterial ecosystem. Most SLOTU of esophageal biota are similar or identical to residents of the upstream oral biota, but the major distinction is that a large majority (82\%) of the esophageal bacteria are known and cultivable. These findings provide evidence for a complex but conserved bacterial population in the normal distal esophagus.
#'
#' (end quote)
#' 
#' A description of the 16S rRNA sequence processing can be found on the mothur-wiki
#' at the link below. A cutoff of 0.10 was used for OTU clustering in that example,
#' and it is taken here as well to create example data, \code{esophagus}, which was 
#' easily imported with the \code{import_mothur()} function.
#'
#' @references 
#' Pei, Z., Bini, E. J., Yang, L., Zhou, M., Francois, F., & Blaser, M. J. (2004). 
#' Bacterial biota in the human distal esophagus.
#' Proceedings of the National Academy of Sciences of the United States of America, 101(12), 4250-4255.
#' \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC384727}
#'
#' mothur-processed files and the sequence data can be downloaded from a zip-file,
#' along with additional description, from the following URL:
#' \url{http://www.mothur.org/wiki/Esophageal_community_analysis}
#' 
#' @name data-esophagus
#' @aliases esophagus
#' @docType data
#' @author Pei et al. \email{zhiheng.pei@@med.nyu.edu}
#' @keywords data
#' @examples
#' ## # Example using esophagus-data in a UniFrac calculation. 
#' ## data(esophagus)
#' ## UniFrac(esophagus, weighted=TRUE)
#' ## UniFrac(esophagus, weighted=FALSE)
#' ## unifrac(t(as(otuTable(esophagus), "matrix")), tre(esophagus) )
#' # # Example importing the mothur example files to create esophagus.
#' # show_mothur_list_cutoffs("~/Dropbox/R/esophagus_example/esophagus.fn.list")
#' # mothlist  <- "~/esophagus_example/esophagus.fn.list"
#' ### mothgroup <- "~/esophagus_example/esophagus.groups"
#' # mothgroup <- "~/esophagus_example/esophagus.good.groups"
#' # mothtree  <- "~/esophagus_example/esophagus.tree"
#' # cutoff    <- "0.10"
#' # esophagus <- import_mothur(mothlist, mothgroup, mothtree, cutoff)
################################################################################
NA
################################################################################
#' (Data) Enterotypes of the human gut microbiome (2011)
#'
#' Published in Nature in early 2011, this work compared (among other things),
#' the faecal microbial communities from 22
#' subjects using complete shotgun DNA sequencing. 
#' Authors further compared these microbial communities with the faecal 
#' communities of subjects from other studies. A total of 280 faecal samples / subjects
#' are represented in this dataset, and 553 genera. The authors claim that the
#' data naturally clumps into three community-level clusters, or ``enterotypes'', 
#' that are not immediately explained by sequencing technology or demographic 
#' features of the subjects, but with potential relevance to understanding 
#' human gut microbiota. 
#'
#' abstract from research article (quoted):
#' 
#' Our knowledge of species and functional composition of the human gut microbiome is rapidly increasing, but it is still based on very few cohorts and little is known about variation across the world. By combining 22 newly sequenced faecal metagenomes of individuals from four countries with previously published data sets, here we identify three robust clusters (referred to as enterotypes hereafter) that are not nation or continent specific. We also confirmed the enterotypes in two published, larger cohorts, indicating that intestinal microbiota variation is generally stratified, not continuous. This indicates further the existence of a limited number of well-balanced host-microbial symbiotic states that might respond differently to diet and drug intake. The enterotypes are mostly driven by species composition, but abundant molecular functions are not necessarily provided by abundant species, highlighting the importance of a functional analysis to understand microbial communities. Although individual host properties such as body mass index, age, or gender cannot explain the observed enterotypes, data-driven marker genes or functional modules can be identified for each of these host properties. For example, twelve genes significantly correlate with age and three functional modules with the body mass index, hinting at a diagnostic potential of microbial markers.
#'
#' (end quote)
#'
#' @references
#' Arumugam, M., Raes, J., Pelletier, E., Le Paslier, D., Yamada, T., Mende, D. R., Fernandes, G. R., et al. (2011).
#' Enterotypes of the human gut microbiome. Nature, 473(7346), 174-180.
#' \url{http://www.nature.com/doifinder/10.1038/nature09944}
#' See supplemental information for subject data. 
#'
#' OTU-clustered data can be downloaded from:
#' http://www.bork.embl.de/Docu/Arumugam_et_al_2011/downloads.html
#'
#' @name data-enterotype
#' @aliases enterotype
#' @docType data
#' @author Arumugam, M., Raes, J., et al.
#' @keywords data
#' @examples
#' # # Try simple network-analysis plot
#' # data(enterotype)
#' # makenetwork(enterotype)
#' # # Filter samples that don't have Enterotype
#' # x <- subset_samples(enterotype, !is.na(Enterotype))
#' # # Create correspondence analysis plot, constrained on the Enterotype category.
#' # calcplot(x ~ Enterotype)
#' # # Alternatively. . .
#' # ent.cca <- cca.phyloseq(x ~ Enterotype)
#' # plot_ordination_phyloseq(ent.cca, x, site_color_category="Enterotype")
#' # # multiple testing of genera correlating with enterotype 2
#' # mt(x, data.frame(sampleData(x))[, "Enterotype"]==2)
#' # # Should return a data.frame, with the following head()
#'                              # # # # # index     teststat   rawp   adjp plower
#' # # # Prevotella                      207 11.469961374 0.0001 0.0088 0.0001
#' # # # Bacteroides                     203 -9.015717540 0.0001 0.0088 0.0001
#' # # # Holdemania                      201 -5.810081084 0.0001 0.0088 0.0001
#' # # # Acetivibrio                     156 -5.246137207 0.0001 0.0088 0.0001
################################################################################
NA
################################################################################
#' (Data) Reproducibility of soil microbiome data (2011)
#'
#' Published in early 2011,
#' this work compared 24 separate soil microbial communities under four treatment
#' conditions via multiplexed/barcoded 454-pyrosequencing of PCR-amplified 16S rRNA gene fragments.
#' The authors found differences in the composition and structure of microbial 
#' communities between soil treatments. 
#' As expected, the soil microbial communities were highly diverse, with a staggering
#' 16,825 different OTUs (species) observed in the included dataset.
#' Interestingly, this study used a larger number of replicates than previous studies of this type,
#' for a total of 56 samples, and the putatively low resampling rate of species 
#' between replicated sequencing trials (``OTU overlap'') was a major concern by
#' the authors.
#'
#' This dataset contains an experiment-level (\code{\link{phyloseq-class}}) object,
#' which in turn contains the taxa-contingency table and soil-treatment table
#' as \code{\link{otuTable-class}} and \code{\link{sampleData-class}} components, respectively.
#'
#' This data was
#' imported from raw files supplied directly by the authors via personal communication
#' for the purposes of including as an example in the \code{\link{phyloseq-package}}. 
#' As this data is sensitive to choices in OTU-clustering parameters, attempts to recreate
#' the \code{otuTable} from the raw sequencing data may give slightly different results
#' than the table provided here. 
#' 
#' abstract from research article (quoted):
#'
#' To determine the reproducibility and quantitation of the amplicon sequencing-based detection approach for analyzing microbial community structure, a total of 24 microbial communities from a long-term global change experimental site were examined. Genomic DNA obtained from each community was used to amplify 16S rRNA genes with two or three barcode tags as technical replicates in the presence of a small quantity (0.1\% wt/wt) of genomic DNA from Shewanella oneidensis MR-1 as the control. The technical reproducibility of the amplicon sequencing-based detection approach is quite low, with an average operational taxonomic unit (OTU) overlap of 17.2\%\code{+/-}2.3\% between two technical replicates, and 8.2\%\code{+/-}2.3\% among three technical replicates, which is most likely due to problems associated with random sampling processes. Such variations in technical replicates could have substantial effects on estimating beta-diversity but less on alpha-diversity. A high variation was also observed in the control across different samples (for example, 66.7-fold for the forward primer), suggesting that the amplicon sequencing-based detection approach could not be quantitative. In addition, various strategies were examined to improve the comparability of amplicon sequencing data, such as increasing biological replicates, and removing singleton sequences and less-representative OTUs across biological replicates. Finally, as expected, various statistical analyses with preprocessed experimental data revealed clear differences in the composition and structure of microbial communities between warming and non-warming, or between clipping and non-clipping. Taken together, these results suggest that amplicon sequencing-based detection is useful in analyzing microbial community structure even though it is not reproducible and quantitative. However, great caution should be taken in experimental design and data interpretation when the amplicon sequencing-based detection approach is used for quantitative analysis of the beta-diversity of microbial communities.
#' 
#' (end quote)
#' 
#' @references Zhou, J., Wu, L., Deng, Y., Zhi, X., Jiang, Y.-H., Tu, Q., Xie, J., et al. 
#'  Reproducibility and quantitation of amplicon sequencing-based detection. 
#'  The ISME Journal. (2011) 5(8):1303-1313. \code{doi:10.1038/ismej.2011.11}
#'
#' The article can be accessed online at \url{http://www.nature.com/ismej/journal/v5/n8/full/ismej201111a.html}
#'
#' @name data-soilrep
#' @aliases soilrep
#' @docType data
#' @author Jizhong Zhou, et al.
#' @keywords data
#' @examples
#' # ################################################################################
#' # # Load the data
#' # ################################################################################
#' # data(soilrep)
#'
#' # ################################################################################
#' # # Richness and sequencing effort example. Accept null hypothesis: 
#' # # No convincing difference in species richness between warmed/unwarmed soils.
#' # ################################################################################
#' # # Build data.frame with total sequencing reads and soil covariates
#' # DF <- data.frame(total.reads=sampleSums(soilrep), as(sampleData(soilrep), "data.frame"))
#' # # Calculate total (estimated) species richness for each sample and combine with DF
#' # DF <- data.frame(DF, t(round(estimateR(t(otuTable(soilrep))))))
#'
#' # # Initialize ggplot data and color layers
#' # mancol <- c(no="blue", yes="red")
#' # p      <- ggplot(DF) + scale_fill_manual(values=mancol) + scale_colour_manual(values=mancol) 
#'
#' # # Build faceted histogram of the total-reads for each treatment type.
#' # p + geom_histogram(aes(x=total.reads, y=..count.., fill=warmed)) + facet_wrap(~ Treatment, 2)
#'
#' # # Plot (estimated) richness versus sequencing effort (total reads)
#' # p + geom_point(aes(x=total.reads, y=S.chao1, color=warmed, shape=clipped), size=3.5)
#'
#' # # Plot (estimated) richness versus observed richness
#' # p + geom_point(aes(x=S.obs, y=S.chao1, color=warmed, shape=clipped), size=3.5)
#'
#' # # Did the warming or clipping treatments affect observed or estimated richness?
#' # par(mfcol=c(2, 2))
#' # boxplot(S.obs ~ warmed, DF, main="Did warming affect observed richness?", xlab="Warmed?")
#' # boxplot(S.chao1 ~ warmed, DF, main="Did warming affect estimated richness (Chao1)?", xlab="Warmed?")
#' # boxplot(S.obs ~ clipped, DF, main="Did clipping affect observed richness?", xlab="Clipped?")
#' # boxplot(S.chao1 ~ clipped, DF, main="Did clipping affect estimated richness (Chao1)?", xlab="Clipped?")
#' # # (For reference, here's a ggplot2 approach to making the boxplot):
#' # p + geom_boxplot(aes(Treatment, S.chao1, color=warmed))
#'
#' # # The treatments do not appear to have affected the
#' # # estimated total richness between warmed/unwarmed soil samples
#' # t.test(x=subset(DF, warmed=="yes")[, "S.chao1"], y=subset(DF, warmed=="no")[, "S.chao1"])
#'
#' # ################################################################################
#' # # A beta diversity comparison.
#' # ################################################################################
#' # jaccdist <- vegdist(t(otuTable(soilrep)), "jaccard")
#'
#' # # soilMDS <- metaMDS(t(otuTable(soilrep)), "jaccard" )
#' # soilMDS <- metaMDS(jaccdist, "jaccard" )
#'
#' # # Add the NMDS coordinates to the soil sample data.frame, DF
#' # DF <- data.frame(DF, scores(soilMDS))
#'
#' # # plot the MDS of jaccard-distances, and shade points by soil treatments
#' # ggplot(DF) + geom_point(aes(x=NMDS1, y=NMDS2, color=Treatment), size=3)
################################################################################
NA
################################################################################