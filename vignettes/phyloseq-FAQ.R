## ---- warning=FALSE, message=FALSE--------------------------------------------
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

## -----------------------------------------------------------------------------
data(esophagus)
plot_tree(esophagus)

## -----------------------------------------------------------------------------
p1 = plot_tree(esophagus, color = "Sample")
p1
p1 + 
  ggtitle("This is my title.") +
  annotate("text", 0.25, 3, 
           color = "orange",
           label = "my annotation")

## -----------------------------------------------------------------------------
data("esophagus")
mdf = psmelt(esophagus)
# Simple bar plot. See plot_bar() for more.
ggplot(mdf, aes(x = Sample, 
                y = Abundance)) + 
  geom_bar(stat = "identity", position = "stack", color = "black")
# Simple heat map. See plot_heatmap() for more.
ggplot(mdf, aes(x = Sample, 
                y = OTU, 
                fill = Abundance)) +
  geom_raster()

