## ----dontrun-basics-vignette, eval=FALSE--------------------------------------
#  vignette("phyloseq-basics")

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library("phyloseq")
library("ggplot2")

## ----ggplot2-themes-----------------------------------------------------------
theme_set(theme_bw())

## -----------------------------------------------------------------------------
data(GlobalPatterns)

## -----------------------------------------------------------------------------
# prune OTUs that are not present in at least one sample
GP <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# Add new human variable to sample data:
sample_data(GP)$human <- factor(human)

## ----richness_estimates0, fig.width=13, fig.height=7--------------------------
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(GP, "human", "SampleType", measures=alpha_meas))

## ----richness_estimates, fig.width=13,height=7--------------------------------
p + geom_boxplot(data=p$data, aes(x=human, y=value, color=NULL), alpha=0.1)

## -----------------------------------------------------------------------------
GP.chl <- subset_taxa(GP, Phylum=="Chlamydiae")

## ----GP-chl-tree, fig.width=15, fig.height=7, message=FALSE, warning=FALSE----
plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="Abundance")

## -----------------------------------------------------------------------------
data(enterotype)

## ----EntAbundPlot, fig.height=6, fig.width=8----------------------------------
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(taxa_sums(enterotype), TRUE)[1:N]/nsamples(enterotype), las=2)

## -----------------------------------------------------------------------------
rank_names(enterotype)

## -----------------------------------------------------------------------------
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10]) 
ent10   <- prune_taxa(TopNOTUs, enterotype)
print(ent10)

## -----------------------------------------------------------------------------
sample_variables(ent10)

## ----entbarplot0, fig.height=6, fig.width=10----------------------------------
plot_bar(ent10, "SeqTech", fill="Enterotype", facet_grid=~Genus)

## ----GPheatmap----------------------------------------------------------------
data("GlobalPatterns")
gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
(p <- plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family"))

## ----GPheatmap-rename-axes----------------------------------------------------
p$scales$scales[[1]]$name <- "My X-Axis"
p$scales$scales[[2]]$name <- "My Y-Axis"
print(p)

## ----plot_sample_network, fig.width=11, fig.height=7, message=FALSE, warning=FALSE----
data(enterotype)
plot_net(enterotype, maxdist=0.4, color="SeqTech", shape="Enterotype")

## ----eval=FALSE---------------------------------------------------------------
#  my.physeq <- import("Biom", BIOMfilename="myBiomFile.biom")
#  my.ord    <- ordinate(my.physeq)
#  plot_ordination(my.physeq, my.ord, color="myFavoriteVarible")

## ----help-import, eval=FALSE--------------------------------------------------
#  help(import)
#  help(ordinate)
#  help(distance)
#  help(plot_ordination)

## ----GP-data-load-------------------------------------------------------------
data(GlobalPatterns)

## ---- eval=FALSE--------------------------------------------------------------
#  GPUF <- UniFrac(GlobalPatterns)

## ----load-precomputed-UF------------------------------------------------------
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))

## -----------------------------------------------------------------------------
GloPa.pcoa = ordinate(GlobalPatterns, method="PCoA", distance=GPUF)

## ----PCoAScree, fig.width=6, fig.height=4-------------------------------------
plot_scree(GloPa.pcoa, "Scree plot for Global Patterns, UniFrac/PCoA")

## ----GPfig5ax1213-------------------------------------------------------------
(p12 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", color="SampleType") + 
  geom_point(size=5) + geom_path() + scale_colour_hue(guide = "none") )
(p13 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", axes=c(1, 3),
  color="SampleType") + geom_line() + geom_point(size=5) )

## ----GP_UF_NMDS0--------------------------------------------------------------
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))
# perform NMDS, set to 2 axes
GP.NMDS <- ordinate(GlobalPatterns, "NMDS", GPUF)
(p <- plot_ordination(GlobalPatterns, GP.NMDS, "samples", color="SampleType") +
  geom_line() + geom_point(size=5) )

## ----GPCAscree0, fig=FALSE----------------------------------------------------
data(GlobalPatterns)
# Take a subset of the GP dataset, top 200 species
topsp <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:200])
GP    <- prune_taxa(topsp, GlobalPatterns)
# Subset further to top 5 phyla, among the top 200 OTUs.
top5ph <- sort(tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
GP     <- subset_taxa(GP, Phylum %in% names(top5ph))
# Re-add human variable to sample data:
sample_data(GP)$human <- factor(human)

## ----GPCAscree, fig.width=8, fig.height=5-------------------------------------
# Now perform a unconstrained correspondence analysis
gpca  <- ordinate(GP, "CCA")
# Scree plot
plot_scree(gpca, "Scree Plot for Global Patterns Correspondence Analysis")

## ----GPCA1234-----------------------------------------------------------------
(p12 <- plot_ordination(GP, gpca, "samples", color="SampleType") + 
  geom_line() + geom_point(size=5) )
(p34 <- plot_ordination(GP, gpca, "samples", axes=c(3, 4), color="SampleType") + 
  geom_line() + geom_point(size=5) )

## ----GPCAspecplot0------------------------------------------------------------
p1  <- plot_ordination(GP, gpca, "species", color="Phylum")
(p1 <- ggplot(p1$data, p1$mapping) + geom_point(size=5, alpha=0.5) + 
  facet_wrap(~Phylum) +  scale_colour_hue(guide = "none") )

## ----GPCAspecplotTopo0--------------------------------------------------------
(p3 <- ggplot(p1$data, p1$mapping) + geom_density2d() +
  facet_wrap(~Phylum) +  scale_colour_hue(guide = "none") )

## ----GPCAjitter0--------------------------------------------------------------
library("reshape2")
# Melt the species-data.frame, DF, to facet each CA axis separately
mdf <- melt(p1$data[, c("CA1", "CA2", "Phylum", "Family", "Genus")], 
            id=c("Phylum", "Family", "Genus") )
# Select some special outliers for labelling
LF <- subset(mdf, variable=="CA2" & value < -1.0)
# build plot: boxplot summaries of each CA-axis, with labels
p <- ggplot(mdf, aes(Phylum, value, color=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(~variable, 2) + 
  scale_colour_hue(guide = "none") +
  theme_bw() + 
  theme( axis.text.x = element_text(angle = -90, vjust = 0.5) )
# Add the text label layer, and render ggplot graphic
(p <- p + geom_text(data=subset(LF, !is.na(Family)),
  mapping = aes(Phylum, value+0.1, color=Phylum, label=Family), 
  vjust=0,
  size=2))

## ----GPtaxaplot0--------------------------------------------------------------
plot_bar(GP, x="human", fill="SampleType", facet_grid= ~ Phylum)

## ----GPdpcoa01----------------------------------------------------------------
# Perform ordination
GP.dpcoa <- ordinate(GP, "DPCoA") 
# Generate default ordination bi-plot
pdpcoa <- 
  plot_ordination(
    physeq = GP, 
    ordination = GP.dpcoa, 
    type="biplot",
    color="SampleType", 
    shape="Phylum")
# Adjust the shape scale manually 
# to make taxa hollow and samples filled (advanced)
shape.fac <- pdpcoa$data$Phylum
man.shapes <- c(19, 21:25)
names(man.shapes) <- c("Samples", levels(shape.fac)[levels(shape.fac)!="Samples"])
p2dpcoa <- pdpcoa + scale_shape_manual(values=man.shapes)
p2dpcoa

## ----GPdpcoa02----------------------------------------------------------------
# Show just Samples or just Taxa
plot_ordination(GP, GP.dpcoa, type="taxa", shape="Phylum")
plot_ordination(GP, GP.dpcoa, type="samples", color="SampleType")
# Split
plot_ordination(GP, GP.dpcoa, type="split",
                color="SampleType", shape="Phylum") +
  ggplot2::scale_colour_discrete()

## ----distancefun--------------------------------------------------------------
data(esophagus)
distance(esophagus, "bray") 
distance(esophagus, "wunifrac") # weighted UniFrac
distance(esophagus, "jaccard") # vegdist jaccard
distance(esophagus, "g") # betadiver method option "g"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  data(esophagus)
#  distance(esophagus, "wUniFrac")
#  distance(esophagus, "uUniFrac")

## -----------------------------------------------------------------------------
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))
# Manually define color-shading vector based on sample type.
colorScale    <- rainbow(length(levels(get_variable(GlobalPatterns, "SampleType"))))
cols          <- colorScale[get_variable(GlobalPatterns, "SampleType")] 
GP.tip.labels <- as(get_variable(GlobalPatterns, "SampleType"), "character")
# This is the actual hierarchical clustering call, specifying average-link clustering
GP.hclust     <- hclust(GPUF, method="average")
plot(GP.hclust, col=cols)

