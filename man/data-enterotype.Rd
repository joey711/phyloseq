\docType{data}
\name{data-enterotype}
\alias{data-enterotype}
\alias{enterotype}
\title{(Data) Enterotypes of the human gut microbiome (2011)}
\description{
  Published in Nature in early 2011, this work compared
  (among other things), the faecal microbial communities
  from 22 subjects using complete shotgun DNA sequencing.
  Authors further compared these microbial communities with
  the faecal communities of subjects from other studies. A
  total of 280 faecal samples / subjects are represented in
  this dataset, and 553 genera. The authors claim that the
  data naturally clumps into three community-level
  clusters, or ``enterotypes'', that are not immediately
  explained by sequencing technology or demographic
  features of the subjects, but with potential relevance to
  understanding human gut microbiota.
}
\details{
  abstract from research article (quoted):

  Our knowledge of species and functional composition of
  the human gut microbiome is rapidly increasing, but it is
  still based on very few cohorts and little is known about
  variation across the world. By combining 22 newly
  sequenced faecal metagenomes of individuals from four
  countries with previously published data sets, here we
  identify three robust clusters (referred to as
  enterotypes hereafter) that are not nation or continent
  specific. We also confirmed the enterotypes in two
  published, larger cohorts, indicating that intestinal
  microbiota variation is generally stratified, not
  continuous. This indicates further the existence of a
  limited number of well-balanced host-microbial symbiotic
  states that might respond differently to diet and drug
  intake. The enterotypes are mostly driven by species
  composition, but abundant molecular functions are not
  necessarily provided by abundant species, highlighting
  the importance of a functional analysis to understand
  microbial communities. Although individual host
  properties such as body mass index, age, or gender cannot
  explain the observed enterotypes, data-driven marker
  genes or functional modules can be identified for each of
  these host properties. For example, twelve genes
  significantly correlate with age and three functional
  modules with the body mass index, hinting at a diagnostic
  potential of microbial markers.

  (end quote)
}
\examples{
# # Try simple network-analysis plot
# data(enterotype)
# makenetwork(enterotype)
# # Filter samples that don't have Enterotype
# x <- subset_samples(enterotype, !is.na(Enterotype))
# # Create correspondence analysis plot, constrained on the Enterotype category.
# calcplot(x ~ Enterotype)
# # Alternatively. . .
# ent.cca <- cca.phyloseq(x ~ Enterotype)
# plot_ordination_phyloseq(ent.cca, x, site_color_category="Enterotype")
# # multiple testing of genera correlating with enterotype 2
# mt(x, data.frame(sampleData(x))[, "Enterotype"]==2)
# # Should return a data.frame, with the following head()
                             # # # # # index     teststat   rawp   adjp plower
# # # Prevotella                      207 11.469961374 0.0001 0.0088 0.0001
# # # Bacteroides                     203 -9.015717540 0.0001 0.0088 0.0001
# # # Holdemania                      201 -5.810081084 0.0001 0.0088 0.0001
# # # Acetivibrio                     156 -5.246137207 0.0001 0.0088 0.0001
}
\author{
  Arumugam, M., Raes, J., et al.
}
\references{
  Arumugam, M., Raes, J., Pelletier, E., Le Paslier, D.,
  Yamada, T., Mende, D. R., Fernandes, G. R., et al.
  (2011). Enterotypes of the human gut microbiome. Nature,
  473(7346), 174-180.
  \url{http://www.nature.com/doifinder/10.1038/nature09944}
  See supplemental information for subject data.

  OTU-clustered data can be downloaded from:
  http://www.bork.embl.de/Docu/Arumugam_et_al_2011/downloads.html
}
\keyword{data}
\keyword{datasets}

