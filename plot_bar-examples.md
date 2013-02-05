
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

plot_bar function: Powerful, flexible phyloseq bar plots
========================================================
The following are examples to help get you started using the `plot_bar` function on your own phyloseq data.

## Global Patterns dataset examples

Load the dataset, and trim to just the *Chlamydiae* phylum.


```r
library("phyloseq")
```


For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.3.12"
```

```r
data("GlobalPatterns")
gp.ch = subset_species(GlobalPatterns, Phylum == "Chlamydiae")
```


### Some Initial Basic Plots
The following is the default barplot when no parameters are given. The dataset is plotted with every sample mapped individually to the horizontal (`x`) axis, and abundance values mapped to the veritcal (`y`) axis. At each sample's horizontal position, the abundance values for each OTU are stacked in order from greatest to least, separate by a thin horizontal line. As long as the parameters you choose to separate the data result in more than one OTU abundance value at the respective position in the plot, the values will be stacked in order as a means of displaying both the sum total value while still representing the individual OTU abundances.


```r
plot_bar(gp.ch)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


Add fill color to represent the Genus to which each OTU belongs.


```r
plot_bar(gp.ch, fill = "Genus")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


Now keep the same fill color, and group the samples together by the `SampleType` variable; essentially, the environment from which the sample was taken and sequenced. 

```r
plot_bar(gp.ch, x = "SampleType", fill = "Genus")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

Note that abundance values for the same OTU from the same `SampleType` will be stacked as separate bar segments, and so the segment lines may not accurately portray the observed richness (because the same OTU might be shown more than once for the same horizontal axis grouping). However, all other aspects of the representation are quantitative, with the total stacked bar height at each horizontal position indicating the sum of all reads for that sample(s). There is not attempt by `plot_bar` to normalize or standardize your data, which is your job to do (using other tools in the phyloseq pacakge, for instance) before attempting to interpret/compare these values between samples.

### More Sophisticated Organization using Facets
In the following example we elected to further organize the data using "facets"  -- separate, adjacent sub-plots. In this case the facets allow us to according to the genus of each OTU. Within each genus facet, the data is further separated by sequencing technology, and the enterotype label for the sample from which each OTU originated is indicated by fill color.


```r
plot_bar(gp.ch, "Family", fill = "Genus", facet_grid = ~SampleType)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



### Further customization using ggplot2 layers

Note that additional customizations of the plot are always possible using standard ggplot2 layers. For example, the following code chunk shows a plot with jittered points add using a second plot layer. 

```r
library("ggplot2")
p = plot_bar(gp.ch, "Family", fill = "Genus", facet_grid = ~SampleType)
p + geom_point(aes(x = Family, y = Abundance), color = "black", position = "jitter", 
    size = 3)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 




## Enterotypes dataset examples

First, load package (if you haven't already), then trim Enterotype data to most abundant 10 genera.


```r
library("phyloseq")
data("enterotype")
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10])
ent10 <- prune_species(TopNOTUs, enterotype)
```


The parameters to `plot_bar` in the following code-chunk were chosen after various trials. We suggest that you also try different parameter settings while you're exploring different features of the data. In addition to the variables names of `sample_data`, the `plot_bar` function recognizes the names of taxonomic ranks, if present. See the help documentation and further details in the examples and on the wiki page. In this example we have also elected to organize data by "facets" (separate, adjacent sub-plots) according to the genus of each OTU. Within each genus facet, the data is further separated by sequencing technology, and the enterotype label for the sample from which each OTU originated is indicated by fill color. Abundance values from different samples and OTUs but having the same variables mapped to the horizontal (`x`) axis are sorted and stacked, with thin horizontal lines designating the boundaries. With this display it is very clear that the choice of sequencing technology had a large effect on which genera were detected, as well as the fraction of OTUs that were assigned to a Genus.


```r
plot_bar(ent10, "SeqTech", fill = "Enterotype", facet_grid = ~Genus)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


You could nix the approach in which OTU abundance values from different samples, different enterotypes, are stacked together and simply shaded differently, and instead opt to separate both the enterotype designation of the samples and the genus designation of the OTUs into one grid. Only a slight modification to the previous function call is necessary in that case (with an added fill to make it even easier to read):


```r
plot_bar(ent10, "Genus", fill = "Genus", facet_grid = SeqTech ~ Enterotype)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 



### Add ggplot2 layer to remove the OTU separation lines

The following example uses more ggplot2-package commands directly for customization, so you need to load the package first (which we did earlier, but I will show it again here for modularity).


```r
library("ggplot2")
```


Now you can save the previous plot as a variable, let's call it `p`, and then add additional ggplot2 layering instructions that will, in effect, remove the dividing lines that separate OTUs from one another in the previous plot. 


```r
p = plot_bar(ent10, "Genus", fill = "Genus", facet_grid = SeqTech ~ Enterotype)
p + geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 

			

## Other key graphics functions in the phyloseq package:

### [plot_ordination](http://joey711.github.com/phyloseq/plot_ordination-examples)

### [plot_heatmap](http://joey711.github.com/phyloseq/plot_heatmap-examples)

### [plot_network](http://joey711.github.com/phyloseq/plot_network-examples)

### [plot_tree](http://joey711.github.com/phyloseq/plot_tree-examples)

### [plot_bar](http://joey711.github.com/phyloseq/plot_bar-examples)

### [plot_richness](http://joey711.github.com/phyloseq/plot_richness-examples)
