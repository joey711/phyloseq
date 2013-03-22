
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

# plot_tree function -
# Powerful tree graphics with ggplot2

This page demos already-constructed examples of phylogenetic trees created via the `plot_tree` function in [the phyloseq package](http://joey711.github.com/phyloseq/), which in turn uses the powerful graphics package called [ggplot2](http://docs.ggplot2.org/current/).

Load the package and datasets


```r
library("phyloseq")
data("esophagus")
data("GlobalPatterns")
```


For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.3.23"
```


We want to plot trees, sometimes even bootstrap values, but notice that the node labels in the `GlobalPatterns` dataset are actually a bit strange. They look like they might be bootstrap values, but they sometimes have two decimals.

```r
head(phy_tree(GlobalPatterns)$node.label, 10)
```

```
##  [1] ""          "0.858.4"   "1.000.154" "0.764.3"   "0.995.2"  
##  [6] "1.000.2"   "0.943.7"   "0.971.6"   "0.766"     "0.611"
```


Could systematically remove the second decimal, but why not just take the first 4 characters instead?

```r
phy_tree(GlobalPatterns)$node.label = substr(phy_tree(GlobalPatterns)$node.label, 
    1, 4)
```


Great, now that we're more happy with the node labels at least looking like bootstrap values, we can move on to using these along with other information about data mapped onto the tree graphic.

The `GlobalPatterns` dataset has many OTUs, more than we would want to try to fit on a tree graphic

```r
ntaxa(GlobalPatterns)
```

```
## [1] 19216
```

So, let's arbitrarily prune to just the first 50 OTUs in `GlobalPatterns`, and store this as `physeq`, which also happens to be the name for most main data parameters of function in the phyloseq package.


```r
physeq = prune_taxa(taxa_names(GlobalPatterns)[1:50], GlobalPatterns)
```


Now let's look at what happens with the default `plot_tree` settings.

```r
plot_tree(physeq)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

By default, black dots are annotated next to tips (OTUs) in the tree, one for each sample in which that OTU was observed. Some have more dots than others. Also by default, the node labels that were stored in the tree were added next to each node without any processing (although we had trimmed their length to 4 characters in the previous step).

What if we want to just see the tree with no sample points next to the tips?

```r
plot_tree(physeq, "treeonly")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

And what about without the node labels either?

```r
plot_tree(physeq, "treeonly", nodeplotblank)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

We can adjust the way branches are rotated to make it look nicer using the `ladderize` parameter.

```r
plot_tree(physeq, "treeonly", nodeplotblank, ladderize = "left")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

```r
plot_tree(physeq, "treeonly", nodeplotblank, ladderize = TRUE)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 

And what if we want to add the OTU labels next to each tip?

```r
plot_tree(physeq, "treeonly", nodeplotblank, label.tips = "taxa_names", ladderize = "left")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


Any `method` parameter argument other than `"sampledodge"` (the default) will not add dodged sample points next to the tips.

```r
plot_tree(physeq, "anythingelse")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 


## Mapping Variables in Data
In the default argument to `method`, `"sampledodge"`, a point is added next to each OTU tip in the tree for every sample in which that OTU was observed. We can then map certain aesthetic features of these points to variables in our data.

### Color
Color is one of the most useful aesthetics in tree graphics that are usually pretty complicated. Color can be mapped to either taxonomic ranks or sample covariates. For instance, we can map color to the type of sample collected (environmental location).

```r
plot_tree(physeq, color = "SampleType", ladderize = "left")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


### node labels
One of the most common reasons to label nodes is to add confidence measures, often a bootstrap value, to the nodes of the tree. The following four graphics show different ways of doing (or not doing) this.

```r
plot_tree(physeq, nodelabf = nodeplotblank, color = "SampleType", ladderize = "left")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-141.png) 

```r
plot_tree(physeq, nodelabf = NULL, color = "SampleType", ladderize = "left")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-142.png) 

```r
plot_tree(physeq, nodelabf = nodeplotboot(), color = "SampleType", ladderize = "left")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-143.png) 

```r
plot_tree(physeq, nodelabf = nodeplotboot(80, 0, 3), color = "SampleType", ladderize = "left")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-144.png) 


### tip labels

```r
plot_tree(physeq, nodelabf = nodeplotboot(80, 0, 3), color = "SampleType", label.tips = "taxa_names", 
    ladderize = "left")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



## The esophagus dataset.
A simple dataset containing tree and OTU-table components.

The esophagus is a small and relatively simple dataset by moderns standards. It only contains 3 samples, no sample-data, and a modest quantity of total sequencing per sample that is a relic of an earlier time when resources for this sort of investigation were sparse and sequencing was expensive. Nevertheless, it makes for a nice sized dataset to start displaying trees. (For more details about the dataset and its origin, try entering `?esophagus` into the command-line once you have loaded the phyloseq package)

The default tree without any additional parameters is plot with black points that indicate in how many different samples each OTU was found. In this case, the term "OTU" is used quite loosely (it is a loose term, after all) to mean entries in the taxononmic data you are plotting; and in the specific case of trees, it means the tips, even if you have already agglomerated the data such that each tip is equivalent to a rank of class or phylum. 


```r
plot_tree(esophagus, title = "Default tree.")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 


If for some reason you just want an unadorned tree, the `"treeonly"` method can be selected. This tends to plot much faster than the annotated tree, and is still a ggplot2 object that you *might* be able to add further layers to manually.


```r
plot_tree(esophagus, "treeonly", title="method = \"treeonly\"")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 


Now let's shade tips according to the sample in which a particular OTU was observed.


```r
plot_tree(esophagus, color = "samples")
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


We can also scale the size of tips according to abundance; usually related to number of sequencing reads, but depends on what you have done with the data prior to this step.


```r
plot_tree(esophagus, size = "abundance")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


Both graphical features included at once.


```r
plot_tree(esophagus, size = "abundance", color = "samples")
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 


There is some overlap of the tip points. Let's adjust the base spacing to spread them out a little bit.


```r
plot_tree(esophagus, size = "abundance", color = "samples", base.spacing = 0.03)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 


Good, now what if we wanted to also display the specific numeric value of OTU abundances that occurred more than 3 times in a given sample? For that, `plot_tree` includes the `min.abundance` parameter, set to `Inf` by default to prevent any point labels from being written.


```r
plot_tree(esophagus, size = "abundance", color = "samples", base.spacing = 0.03, 
    min.abundance = 3)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 


# More Examples with the Global Patterns dataset

Subset Global Patterns dataset to just the observed Archaea.


```r
gpa <- subset_species(GlobalPatterns, Kingdom == "Archaea")
```


The number of different Archaeal species from this dataset is small enough for a decent tree plot.


```r
nspecies(gpa)
```

```
## [1] 208
```


That is to say, it is reasonable to consider displaying the phylogenetic tree directly for `gpa`. Too many OTUs means a tree that is pointless to attempt to display in its entirety in one graphic of a standard size and font. So the whole `GlobalPatterns` dataset probably a bad idea.


```r
nspecies(GlobalPatterns)
```

```
## [1] 19216
```


Some patterns are immediately discernable with minimal parameter choices:


```r
plot_tree(gpa, color = "SampleType")
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26.png) 




```r
plot_tree(gpa, color = "Phylum")
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.png) 




```r
plot_tree(gpa, color = "SampleType", shape = "Phylum")
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28.png) 




```r
plot_tree(gpa, color = "Phylum", label.tips = "Genus")
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29.png) 


However, the text-label size scales with number of species, and with common graphics-divice sizes/resolutions, these ~200 taxa still make for a somewhat crowded graphic. 

Let's instead subset further to just the Crenarchaeota


```r
gpac <- subset_species(gpa, Phylum == "Crenarchaeota")
plot_tree(gpac, color = "SampleType", shape = "Genus")
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30.png) 




```r
plot_tree(gpac, color = "SampleType", label.tips = "Genus")
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31.png) 


Let's add some abundance information. Notice that the default spacing gets a little crowded when we map species-abundance to point-size:


```r
plot_tree(gpac, color = "SampleType", shape = "Genus", size = "abundance", plot.margin = 0.4)
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32.png) 


So let's spread it out a little bit with the `base.spacing` parameter, and while we're at it, let's call off the node labels...


```r
plot_tree(gpac, nodelabf = nodeplotblank, color = "SampleType", shape = "Genus", 
    size = "abundance", base.spacing = 0.04, plot.margin = 0.4)
```

![plot of chunk unnamed-chunk-33](figure/unnamed-chunk-33.png) 


## Chlamydiae-only tree


```r
GP.chl <- subset_species(GlobalPatterns, Phylum == "Chlamydiae")
plot_tree(GP.chl, color = "SampleType", shape = "Family", label.tips = "Genus", 
    size = "abundance", plot.margin = 0.6)
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34.png) 

			

## Other key graphics functions in the phyloseq package:

### [plot_ordination](http://joey711.github.com/phyloseq/plot_ordination-examples)

### [plot_heatmap](http://joey711.github.com/phyloseq/plot_heatmap-examples)

### [plot_network](http://joey711.github.com/phyloseq/plot_network-examples)

### [plot_tree](http://joey711.github.com/phyloseq/plot_tree-examples)

### [plot_bar](http://joey711.github.com/phyloseq/plot_bar-examples)

### [plot_richness](http://joey711.github.com/phyloseq/plot_richness-examples)
