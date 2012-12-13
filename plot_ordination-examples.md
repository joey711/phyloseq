
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>


# plot_ordination

## see also...
The operation of this function also depends a lot on the 

[distance](https://github.com/joey711/phyloseq/wiki/distance)

and

[ordinate](https://github.com/joey711/phyloseq/wiki/ordinate)

functions. See their wiki-pages for further details and examples.

## Examples

Load the necessary packages and data.

```r
library("phyloseq")
library("ggplot2")
data(GlobalPatterns)
```


"Trim" data. This is useful for plotting, and in this case, also useful for making examples that run in a short amount of time. Your reasoning and decisions in trimming are extremely important, and up to you. I am using several different methods of trimming here, for illustration and because the extent of data reduction is useful for my purposes. However, I make no assertion that these are the "right" approach(es) for your data, but rather, I highly recommend that you think hard about any trimming you do, and only commit to including it in your final analysis pipeline if you can defend the choices and have checked that they are robust. 

Need to clean the zeros from GlobalPatterns. While we're at it, let's remove any OTUs observed less than 5 times, cumulative for all samples


```r
GP <- prune_taxa(speciesSums(GlobalPatterns) > 5, GlobalPatterns)
nspecies(GP)
```

```
## [1] 13711
```


Still more than $13000$ OTUs. Let's filter taxa that don't show up at least 5 times in 5 or more samples.


```r
wh0 <- genefilterSample(GP, filterfunSample(function(x) {
    x > 5
}), A = 5)
GP <- prune_taxa(wh0, GP)
```


Do an additional trimming by cumulative abundance of phyla, take only the top 5


```r
phylum.sum <- tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE)
top5phyla <- names(sort(phylum.sum, TRUE))[1:5]
GP1 <- subset_taxa(GP, Phylum %in% top5phyla)
```


We will want to investigate a major prior among the samples, which is that some are human-associated microbiomes, and some are not. Define a human-associated versus non-human categorical variable:


```r
human <- getVariable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sampleData(GP1)$human <- factor(human)
```




## Examples using the Global Patterns dataset and `plot_ordination`

Let's start by plotting just the species/taxa, and shading the points by Phylum. Note that even in our "trimmed" dataset there are `nspecies(GP1)`=2265 OTUs.


```r
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



This is a complicated looking plot, but that's not necessarily good. There is actually a lot of overplotting/occlusion, which means that the high number of points is getting in the way of our visual understanding of the data. There are several ways to deal with this in ggplot2, for example, facetting:


```r
p1 + facet_wrap(~Phylum, 3)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



Next, let's plot only the samples, and shade the points by "SampleType" while also modifying the shape according to whether they are human-associated. There are a few additional ggplot2 layers added to make the plot even nicer...


```r
p2 = plot_ordination(GP1, GP.ord, type = "samples", color = "SampleType", shape = "human")
p2 + geom_line() + geom_point(size = 5) + ggtitle("samples")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



Now let's try combining both the samples and species together in one "biplot".


```r
p3 = plot_ordination(GP1, GP.ord, type = "biplot", color = "SampleType", shape = "Phylum", 
    title = "biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names <- getTaxa(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values = GP1.shape)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



Hmmm, the overlap problem is affecting this again, let's try the "split" organization of the biplot, in which the samples/species are separated on two panels...

```r
p4 = plot_ordination(GP1, GP.ord, type = "split", color = "Phylum", shape = "human", 
    label = "SampleType", title = "split")
# Adjust colors to make species black The following function reproduces
# ggplot2's default color scale.  From:
# http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
color.names <- levels(p4$data$Phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
p4 + scale_color_manual(values = p4cols)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 




## Other key graphics functions in the phyloseq package:

### [plot_ordination](http://joey711.github.com/phyloseq/plot_ordination-examples)

### [plot_tree](http://joey711.github.com/phyloseq/plot_tree-examples)

### [plot_network](http://joey711.github.com/phyloseq/plot_network-examples)

### [plot_richness](http://joey711.github.com/phyloseq/plot_richness-examples)

### [plot_bar](http://joey711.github.com/phyloseq/plot_bar-examples)
