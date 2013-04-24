
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# plot_heatmap
A tutorial with Examples


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.5'
```

```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```


ggplot2 package theme set. See [the ggplot2 online documentation](http://docs.ggplot2.org/current/) for further help.


```r
theme_set(theme_bw())
```



The following demonstrates some uses of the `plot_heatmap` function in [the phyloseq package](http://joey711.github.com/phyloseq/) for R and Bioconductor.

## Introduction - Create an ecologically-organized heatmap 
 
In a [2010 article in BMC Genomics](http://www.biomedcentral.com/1471-2105/11/45), Rajaram and Oono describe an approach to creating a heatmap using ordination methods (namely, NMDS and PCA) to organize the rows and columns instead of (hierarchical) cluster analysis. In many cases the ordination-based ordering does a much better job than h-clustering at providing an order of elements that is easily interpretable. The authors provided an immediately useful example of their approach as [the NeatMap package for R](http://cran.r-project.org/web/packages/NeatMap/index.html). The NeatMap package can be used directly on the abundance table (`"otuTable"`-class) of phylogenetic-sequencing data, but the NMDS or PCA ordination options that it supports are not based on ecological distances. To fill this void, and because phyloseq already provides support for a large number of [ecological distances](http://joey711.github.io/phyloseq/distance) and [ordination methods](http://joey711.github.io/phyloseq/ordinate), phyloseq now includes the `plot_heatmap()` function: an ecology-oriented variant of the NeatMap approach to organizing a heatmap and build it using ggplot2 graphics tools. The [distance](http://joey711.github.io/phyloseq/distance) and [method](http://joey711.github.io/phyloseq/ordinate) arguments are the same as for the [plot_ordination](http://joey711.github.io/phyloseq/plot_ordination) function, and support large number of distances and ordination methods, respectively, with a strong leaning toward ecology. This function also provides the options to re-label the OTU and sample axis-ticks with a taxonomic name and/or sample variable, respectively, in the hope that this might hasten your interpretation of the patterns (See the documentation for the `sample.label` and `taxa.label` arguments, and the examples below). Note that this function makes no attempt to overlay dendrograms from hierarchical clustering next to the axes, as hierarchical clustering is not used to organize these plots. Also note that each re-ordered axis repeats at the edge, and so apparent clusters at the far right/left or top/bottom of the heat-map may actually be the same. For now, the placement of this edge can be considered arbitrary, so beware of this artifact of the graphic and visually check if there are two "mergeable" clusters at the edges of a particular axis. If you benefit from this phyloseq-specific implementation of [the NeatMap approach](http://cran.r-project.org/web/packages/NeatMap/index.html), please cite [the NeatMap article](http://www.biomedcentral.com/1471-2105/11/45), as well as phyloseq.

## Heatmap colors don't have to be so hot

Traditionally heatmaps have been used to emphasize data that is above or below a threshold as "hot" or "cold" colors, respectively. In the case of OTU-abundance data, however, it seems the most common need is to see the relative patterns of high-abundance OTUs against a background of taxa that are mostly low-abundance or absent in a sparse matrix. Furthermore, there is usually not an obvious or intrinsically meaningful abundance value to use as a suitable threshold for the traditional "cold/hot" display. For these reasons, the default color scheme in `plot_heatmap` maps a very dark blue color to the lowest abundance values, up to a very light blue for the highest abundance values. The dark blue for the lowest abundance values is not very much lighter than black - the color used to represent missing or zero abundance values by default - for a coherent, blue-oriented color scheme in which the eye should be drawn to the lighter shades.

If, for whatever reason, you need to change this default color scheme, it is possible through the `low`, `high`, and `na.value` arguments. Several examples are provided below. The character-string values supplied to these arguments need to be the names of R colors. There are over 600 English color names that are understood by R (try `colors()` at the R terminal), as well as other finely-resolved color gradient nomenclatures. The examples below use a 6-digit hexadecimal color representation, and a nice [table summary of these colors is available at the R Cookbook](http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2\)/).

For further details, not that the `high`, `low`, and `na.value` parameters are passed along to [ggplot2's scale_gradient function](http://had.co.nz/ggplot2/scale_gradient.html), which does an excellent job selecting suitable colors of your gradient, provided that you select colors that make sense to have at two ends of a gradient.

I also got some useful ideas and suggestions at the following [WordPress page regarding the construction of heatmaps using ggplot2's](http://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/) `geom_tile`.

## Adjust color scale of your `plot_heatmap`

By default, the color mapping of `plot_heatmap` is transformed to a log-scale of base 4, using `log_trans(4)` from [the scales package](http://cran.r-project.org/web/packages/scales/index.html). This is an arbitrary choice that you might need to adjust based on your needs and data. If specifying an alternative transformation object to the `trans` argument, you probably need to load the scales package first. Since [scales](http://cran.r-project.org/web/packages/scales/index.html) is a required package for `phyloseq`, you should already have it installed if you are at this point. Any transformation object that is valid for the scales package should work here, but the relative contrast and the way it represents your data could change dramatically based on this choice, so make this selection carefully; or better yet, try several different transformations if you think data is being "left in the background" or too much information is being "pushed to the foreground", for example.

## Plot a 300-taxa dataset
The following two lines subset the dataset to just the top 300 most abundant Bacteria taxa across all samples (in this case, with no prior preprocessing. Not recommended, but quick).


```r
gpt <- subset_taxa(GlobalPatterns, Kingdom == "Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt), TRUE)[1:300]), gpt)
plot_heatmap(gpt, sample.label = "SampleType")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


## Subset a smaller dataset based on an Archaeal phylum

Subset the dataset to something manageable that can be reasonably
represented in one plot. In the following examples, the Crenarchaeota phylum.


```r
gpac <- subset_taxa(GlobalPatterns, Phylum == "Crenarchaeota")
```



## Default `plot_heatmap` settings

Now let's see how our `plot_heatmap` function works with all default settings.

```r
plot_heatmap(gpac)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 



## Re-label by a sample variable and taxonomic family

Here is an example re-labelling based on the "SampleType" sample variable and the taxonomic rank of "Family".


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 



### Now repeat the plot, but change the color scheme.

Changing the color scheme might be worthwhile, depending on the graphics device or paper
on which you want to display the heatmap. 


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low = "#000033", 
    high = "#CCFF66")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



Here is a dark-blue to red scheme.


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low = "#000033", 
    high = "#FF3300")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



A very dark-blue to very light-blue scheme


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low = "#000033", 
    high = "#66CCFF")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



Here is a "dark on light" color scheme. Note that we change the background value
(the value of the NA and 0 elements)


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low = "#66CCFF", 
    high = "#000033", na.value = "white")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



This is a similar color scheme as the previous, but the "near zero" color is closer to 
a cream color, and the colors themselves are closer to blue-grey. This is better
overall contrast than a lot of schemes, but may not be as exciting.


```r
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low = "#FFFFCC", 
    high = "#000033", na.value = "white")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 



## Now try different ordination methods, distances

Now try the default color scheme, but using different ecological distances/ordinations.
For example, NMDS ordination on the jaccard distance.


```r
plot_heatmap(gpac, "NMDS", "jaccard")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 



Detrended correspondence analysis.


```r
plot_heatmap(gpac, "DCA", "none", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 



Unconstrained redundancy analysis (Principle Components Analysis, PCA)


```r
plot_heatmap(gpac, "RDA", "none", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 



PCoA/MDS ordination on the (default) bray-curtis distance.


```r
plot_heatmap(gpac, "PCoA", "bray", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 



MDS/PCoA ordination on the Unweighted-UniFrac distance.


```r
plot_heatmap(gpac, "PCoA", "unifrac", "SampleType", "Family")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



Now try weighted-UniFrac distance and MDS/PCoA ordination.


```r
plot_heatmap(gpac, "MDS", "unifrac", "SampleType", "Family", weighted = TRUE)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 



Here is how you might create a heatmap using base-R graphics
and the more common (but problematic) hierarchical clustering
organization, in case you want to compare with `plot_heatmap`, for example.


```r
heatmap(otuTable(gpac))
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 

			


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


