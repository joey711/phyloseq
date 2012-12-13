
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

# plot_tree function -
# Powerful tree graphics with ggplot2

This page contains many already-constructed examples of trees created by the `plot_tree` function in the phyloseq package, which in turn is uses the powerful graphics framework package called [ggplot2](http://docs.ggplot2.org/current/).

## The esophagus dataset.
A simple dataset containing tree and OTU-table components.

The esophagus is a small and relatively simple dataset by moderns standards. It only contains 3 samples, no sample-data, and a modest quantity of total sequencing per sample that is a relic of an earlier time when resources for this sort of investigation were sparse and sequencing was expensive. Nevertheless, it makes for a nice sized dataset to start displaying trees. (For more details about the dataset and its origin, try entering `?esophagus` into the command-line once you have loaded the phyloseq package)

Load the package and dataset


```r
library("phyloseq")
data("esophagus")
```


The default tree without any additional parameters is plot with black points that indicate in how many different samples each OTU was found. In this case, the term "OTU" is used quite loosely (it is a loose term, after all) to mean entries in the taxononmic data you are plotting; and in the specific case of trees, it means the tips, even if you have already agglomerated the data such that each tip is equivalent to a rank of class or phylum. 


```r
plot_tree(esophagus, title = "Default tree.")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


If for some reason you just want a completely unadorned tree, the `"treeonly"` method can be selected. This tends to plot much faster than the annotated tree, and is still a ggplot2 object that you *might* be able to add further layers to manually.


```r
plot_tree(esophagus, "treeonly", title="method = \"treeonly\"")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


Now let's shade tips according to the sample in which a particular OTU was observed.


```r
plot_tree(esophagus, color = "samples")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


We can also scale the size of tips according to abundance; usually related to number of sequencing reads, but depends on what you have done with the data prior to this step.


```r
plot_tree(esophagus, size = "abundance")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Both graphical features included at once.


```r
plot_tree(esophagus, size = "abundance", color = "samples")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


There is some overlap of the tip points. Let's adjust the base spacing to spread them out a little bit.


```r
plot_tree(esophagus, size = "abundance", color = "samples", base.spacing = 0.03)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Good, now what if we wanted to also display the specific numeric value of OTU abundances that occurred more than 3 times in a given sample? For that, `plot_tree` includes the `min.abundance` parameter, set to `Inf` by default to prevent any point labels from being written.


```r
plot_tree(esophagus, size = "abundance", color = "samples", base.spacing = 0.03, 
    min.abundance = 3)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


# Examples with the Global Patterns dataset

Load the dataset (quotes are optional, did you know that?)


```r
data(GlobalPatterns)
```


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

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 




```r
plot_tree(gpa, color = "Phylum")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 




```r
plot_tree(gpa, color = "SampleType", shape = "Phylum")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 




```r
plot_tree(gpa, color = "Phylum", label.tips = "Genus")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 


However, the text-label size scales with number of species, and with common graphics-divice sizes/resolutions, these ~200 taxa still make for a somewhat crowded graphic. 

Let's instead subset further to just the Crenarchaeota


```r
gpac <- subset_species(gpa, Phylum == "Crenarchaeota")
plot_tree(gpac, color = "SampleType", shape = "Genus")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 




```r
plot_tree(gpac, color = "SampleType", label.tips = "Genus")
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


Let's add some abundance information. Notice that the default spacing gets a little crowded when we map species-abundance to point-size:


```r
plot_tree(gpac, color = "SampleType", shape = "Genus", size = "abundance")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


So let's spread it out a little bit with the base.spacing parameter.


```r
plot_tree(gpac, color = "SampleType", shape = "Genus", size = "abundance", base.spacing = 0.05)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 


## Chlamydiae-only tree


```r
GP.chl <- subset_species(GlobalPatterns, Phylum == "Chlamydiae")
plot_tree(GP.chl, color = "SampleType", shape = "Family", label.tips = "Genus", 
    size = "abundance")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 

			

## Other key graphics functions in the phyloseq package:

### [plot_ordination](http://joey711.github.com/phyloseq/plot_ordination-examples)

### [plot_heatmap](http://joey711.github.com/phyloseq/plot_heatmap-examples)

### [plot_network](http://joey711.github.com/phyloseq/plot_network-examples)

### [plot_tree](http://joey711.github.com/phyloseq/plot_tree-examples)

### [plot_bar](http://joey711.github.com/phyloseq/plot_bar-examples)

### [plot_richness](http://joey711.github.com/phyloseq/plot_richness-examples)
