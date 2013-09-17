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
#' # Example using esophagus-data in a UniFrac calculation. 
#' data(esophagus)
#' UniFrac(esophagus, weighted=TRUE)
#' UniFrac(esophagus, weighted=FALSE)
#' #
#' # How to re-create the esophagus dataset using import_mothur function
#' mothlist  <- system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
#' mothgroup <- system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
#' mothtree  <- system.file("extdata", "esophagus.tree.gz", package="phyloseq")
#' show_mothur_list_cutoffs(mothlist)
#' cutoff    <- "0.10"
#' esophman  <- import_mothur(mothlist, mothgroup, mothtree, cutoff)	
#' identical(esophagus, esophman)
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
#' Arumugam, M., et al. (2011). Enterotypes of the human gut microbiome.
#'
#' Nature, 473(7346), 174-180.
#'
#' \url{http://www.nature.com/doifinder/10.1038/nature09944}
#' See supplemental information for subject data. 
#'
#' OTU-clustered data was downloaded from the publicly-accessible:
#'
#' \url{http://www.bork.embl.de/Docu/Arumugam_et_al_2011/downloads.html}
#'
#' @name data-enterotype
#' @aliases enterotype
#' @docType data
#' @author Arumugam, M., Raes, J., et al.
#' @keywords data
#' @examples
#' # Try simple network-analysis plot
#' data(enterotype)
#' ig <- make_network(enterotype, "samples", max.dist=0.3)
#' plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#' # 
#' # Filter samples that don't have Enterotype
#' x <- subset_samples(enterotype, !is.na(Enterotype))
#' # 
#' # Alternatively. . .
#' ent.cca <- ordinate(x ~ Enterotype, "CCA")
#' plot_ordination(x, ent.cca, color="Enterotype")
#' plot_ordination(x, ent.cca, "biplot")
#' plot_ordination(x, ent.cca, "split", color="Enterotype")
#' # 
#' # # multiple testing of genera correlating with enterotype 2
#' # mt(x, data.frame(sample_data(x))[, "Enterotype"]==2)
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
#' as \code{\link{otu_table-class}} and \code{\link{sample_data-class}} components, respectively.
#'
#' This data was
#' imported from raw files supplied directly by the authors via personal communication
#' for the purposes of including as an example in the \code{\link{phyloseq-package}}. 
#' As this data is sensitive to choices in OTU-clustering parameters, attempts to recreate
#' the \code{otu_table} from the raw sequencing data may give slightly different results
#' than the table provided here. 
#' 
#' abstract from research article (quoted):
#'
#' To determine the reproducibility and quantitation of the amplicon sequencing-based 
#' detection approach for analyzing microbial community structure, a total of 24 microbial
#' communities from a long-term global change experimental site were examined. Genomic DNA
#' obtained from each community was used to amplify 16S rRNA genes with two or three 
#' barcode tags as technical replicates in the presence of a small quantity (0.1\% wt/wt)
#' of genomic DNA from Shewanella oneidensis MR-1 as the control. The technical 
#' reproducibility of the amplicon sequencing-based detection approach is quite low, 
#' with an average operational taxonomic unit (OTU) overlap of 17.2\%\code{+/-}2.3\% 
#' between two technical replicates, and 8.2\%\code{+/-}2.3\% among three technical 
#' replicates, which is most likely due to problems associated with random sampling processes.
#' Such variations in technical replicates could have substantial effects on estimating 
#' beta-diversity but less on alpha-diversity. A high variation was also observed in the 
#' control across different samples (for example, 66.7-fold for the forward primer), 
#' suggesting that the amplicon sequencing-based detection approach could not be quantitative.
#' In addition, various strategies were examined to improve the comparability of amplicon 
#' sequencing data, such as increasing biological replicates, and removing singleton sequences 
#' and less-representative OTUs across biological replicates. Finally, as expected, various 
#' statistical analyses with preprocessed experimental data revealed clear differences in 
#' the composition and structure of microbial communities between warming and non-warming, 
#' or between clipping and non-clipping. Taken together, these results suggest that amplicon 
#' sequencing-based detection is useful in analyzing microbial community structure even 
#' though it is not reproducible and quantitative. However, great caution should be taken 
#' in experimental design and data interpretation when the amplicon sequencing-based detection 
#' approach is used for quantitative analysis of the beta-diversity of microbial communities.
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
#' # Load the data
#' data(soilrep)
#' ################################################################################
#' # Alpha diversity (richness) example. Accept null hypothesis: 
#' # No convincing difference in species richness between warmed/unwarmed soils.
#' ################################################################################
#' # Graphically compare richness between the different treatments.
#' man.col <- c(WC="red", WU="brown", UC="blue", UU="darkgreen")
#' (p <- plot_richness_estimates(soilrep, x="Treatment", color="Treatment") )
#' # Add boxplots using ggplot2
#' # library(ggplot2)
#' # p + geom_boxplot() + scale_color_manual(values=man.col)
#' # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#' # The treatments do not appear to have affected the
#' # estimated total richness between warmed/unwarmed soil samples
#' # Test this formally:
#' DF <- data.frame(sample_data(soilrep), estimate_richness(soilrep) )
#' t.test(x=subset(DF, warmed=="yes")[, "Chao1"], y=subset(DF, warmed=="no")[, "Chao1"])
#' ################################################################################
#' # A beta diversity comparison.
#' ################################################################################
#' # Perform non-metric multidimensional scaling, using Bray-Curtis distance
#' soil.NMDS <- ordinate(soilrep, "NMDS", "bray")
#' (p <- plot_ordination(soilrep, soil.NMDS, "samples", color="Treatment") )
#' # Additional formatting using ggplot2
#' # library(ggplot2)
#' # ( p <- p + geom_point(size=5, alpha=0.5) + facet_grid(warmed ~ clipped) )
################################################################################
NA
################################################################################
################################################################################
#' (Data) Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample (2011)
#'
#' Published in PNAS in early 2011. This work compared the microbial 
#' communities from 25 environmental samples and three known ``mock communities''
#' -- a total of 9 sample types -- at a depth averaging 3.1 million reads per sample.
#' Authors were able to reproduce diversity patterns seen in many other 
#' published studies, while also invesitigating technical issues/bias by 
#' applying the same techniques to simulated microbial communities of known
#' composition.
#'
#' abstract from research article (quoted):
#' 
#' The ongoing revolution in high-throughput sequencing continues to democratize the ability of small groups of investigators to map the microbial component of the biosphere. In particular, the coevolution of new sequencing platforms and new software tools allows data acquisition and analysis on an unprecedented scale. Here we report the next stage in this coevolutionary arms race, using the Illumina GAIIx platform to sequence a diverse array of 25 environmental samples and three known ``mock communities'' at a depth averaging 3.1 million reads per sample. We demonstrate excellent consistency in taxonomic recovery and recapture diversity patterns that were previously reported on the basis of metaanalysis of many studies from the literature (notably, the saline/nonsaline split in environmental samples and the split between host-associated and free-living communities). We also demonstrate that 2,000 Illumina single-end reads are sufficient to recapture the same relationships among samples that we observe with the full dataset. The results thus open up the possibility of conducting large-scale studies analyzing thousands of samples simultaneously to survey microbial communities at an unprecedented spatial and temporal resolution.
#'
#' (end quote)
#'
#' Many thanks to J. Gregory Caporaso for directly providing the OTU-clustered data files
#' for inclusion in this package.
#'
#' @references
#' Caporaso, J. G., et al. (2011). 
#' Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample.
#' PNAS, 108, 4516-4522.
#' PMCID: PMC3063599
#'
#' The primary article can be viewed/downloaded at:
#' \url{http://www.pnas.org/content/108/suppl.1/4516.short}
#'
#' @name data-GlobalPatterns
#' @aliases GlobalPatterns
#' @docType data
#' @author Caporaso, J. G., et al.
#' @keywords data
#'
#' @seealso 
#'  The examples on the phyloseq wiki page for \code{\link{plot_ordination}} show 
#'  many more examples:
#'
#' \url{https://github.com/joey711/phyloseq/wiki/plot_ordination}
#'
#' @examples
#' data(GlobalPatterns)
#' # Remove unobserved taxa
#' GP0   <- prune_species(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
#' # Perform ordination (in this case, detrended correspondence analysis)
#' gpdca <- ordinate(GP0, "DCA")
#' # Create plot of samples
#' plot_ordination(GP0, gpdca, color="SampleType", title="DCA on abundances, first two axes")
#' # # More complicated plot facetting by phylum.
#' # library("ggplot2")
#' # plot_ordination(GP0, gpdca, color="SampleType", title="DCA on abundances, first two axes") + geom_line()
#' # plot_ordination(GP0, gpdca, "taxa", color="Kingdom") + facet_wrap(~Phylum, 8)
################################################################################
NA
################################################################################
