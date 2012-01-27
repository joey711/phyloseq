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