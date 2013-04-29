
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# plot_ordination examples

---

The operation of the `plot_ordination` function also depends a lot on the 

### [distance](http://joey711.github.io/phyloseq/distance)

and

### [ordinate](http://joey711.github.io/phyloseq/ordinate)

functions. See their tutorials for further details and examples.

Also, the phyloseq package includes a "convenience function" for subsetting from large collections of points in an ordination, called `subset_ord_plot`. There is a separate [subset_ord_plot tutorial](http://joey711.github.com/phyloseq/subset_ord_plot-examples) for further details and examples. 


## Examples

Load the necessary packages and data.

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.6'
```

```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```

```r
data(GlobalPatterns)
```


ggplot2 package theme set. See [the ggplot2 online documentation](http://docs.ggplot2.org/current/) for further help.


```r
theme_set(theme_bw())
```


"Trim" data. This is useful for plotting, and in this case, also useful for making examples that run in a short amount of time. Your reasoning and decisions in trimming are extremely important, and up to you. I am using several different methods of trimming here, for illustration and because the extent of data reduction is useful for my purposes. However, I make no assertion that these are the "right" approach(es) for your data, but rather, I highly recommend that you think hard about any trimming you do, and only commit to including it in your final analysis pipeline if you can defend the choices and have checked that they are robust. 

Need to clean the zeros from GlobalPatterns. While we're at it, let's remove any OTUs observed less than 5 times, cumulative for all samples


```r
GP <- prune_taxa(taxa_sums(GlobalPatterns) > 5, GlobalPatterns)
ntaxa(GP)
```

```
## [1] 13711
```


Still more than 13711 OTUs. Let's filter taxa that don't show up at least 5 times in 5 or more samples.


```r
wh0 <- genefilter_sample(GP, filterfunSample(function(x) {
    x > 5
}), A = 5)
GP <- prune_taxa(wh0, GP)
```


Do an additional trimming by cumulative abundance of phyla, take only the top 5


```r
phylum.sum <- tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum, na.rm = TRUE)
top5phyla <- names(sort(phylum.sum, TRUE))[1:5]
GP1 <- subset_taxa(GP, Phylum %in% top5phyla)
```


We will want to investigate a major prior among the samples, which is that some are human-associated microbiomes, and some are not. Define a human-associated versus non-human categorical variable:


```r
human <- get_variable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP1)$human <- factor(human)
```




## Examples using the Global Patterns dataset and `plot_ordination`

### Just OTUs

Let's start by plotting just the OTUs, and shading the points by Phylum. Note that even in our "trimmed" dataset there are `ntaxa(GP1)=` 2265 OTUs.


```r
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
```

```
## Warning: is.na() applied to non-(list or vector) of type 'NULL'
```

```r
print(p1)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



This is a complicated looking plot, but that's not necessarily good. There is actually a lot of overplotting/occlusion, which means that the high number of points is getting in the way of our visual understanding of the data. There are several ways to deal with this in ggplot2, for example, facetting:


```r
p1 + facet_wrap(~Phylum, 3)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


### Just samples

Next, let's plot only the samples, and shade the points by "SampleType" while also modifying the shape according to whether they are human-associated. There are a few additional ggplot2 layers added to make the plot even nicer...


```r
p2 = plot_ordination(GP1, GP.ord, type = "samples", color = "SampleType", shape = "human")
```

```
## Warning: is.na() applied to non-(list or vector) of type 'NULL'
```

```r
p2 + geom_polygon(aes(fill = SampleType)) + geom_point(size = 5) + ggtitle("samples")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



### Combined/biplot graphic

The `plot_ordination` function can also automatically create two different graphic layouts in which both the samples and OTUs are plotted together in one "biplot". Note that this requires methods that are not intrinsically samples-only ordinations. For example, this doesn't work with UniFrac/PCoA.


```r
p3 = plot_ordination(GP1, GP.ord, type = "biplot", color = "SampleType", shape = "Phylum", 
    title = "biplot")
```

```
## Warning: is.na() applied to non-(list or vector) of type 'NULL'
```

```r
# Some stuff to modify the automatic shape scale
GP1.shape.names <- getTaxa(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values = GP1.shape)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



Hmmm, the overlap problem is affecting this again, let's try the "split" organization of the biplot, in which the samples/OTUs are separated on two panels...

```r
p4 = plot_ordination(GP1, GP.ord, type = "split", color = "Phylum", shape = "human", 
    label = "SampleType", title = "split")
```

```
## Warning: is.na() applied to non-(list or vector) of type 'NULL'
```

```r
# Adjust colors to make OTUs black The following function reproduces
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

			


---

### Other tutorial pages for the phyloseq package:

#### [distance](distance.html)

#### [Example-Data](Example-Data.html)

#### [future-devel](future-devel.html)

#### [import-data](import-data.html)

#### [index](index.html)

#### [install](install.html)

#### [merge](merge.html)

#### [plot_bar-examples](plot_bar-examples.html)

#### [plot_heatmap-examples](plot_heatmap-examples.html)

#### [plot_network-examples](plot_network-examples.html)

#### [plot_ordination-examples](plot_ordination-examples.html)

#### [plot_richness-examples](plot_richness-examples.html)

#### [plot_tree-examples](plot_tree-examples.html)

#### [preprocess](preprocess.html)

#### [rebuild-all-from-Rmarkdown.R](rebuild-all-html-from-Rmarkdown.R)

#### [subset_ord_plot-examples](subset_ord_plot-examples.html)

#### [tutorials-index](tutorials-index.html)


