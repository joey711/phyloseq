
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

plot_network Examples
========================================================
## Load phyloseq and the "enterotypes" dataset


```r
library(phyloseq)
library(ggplot2)
data(enterotype)
```

For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.5.4"
```

```r
packageDescription("ggplot2")$Version
```

```
## [1] "0.9.3.1"
```


Create an igraph-based network based on the default distance method, "Jaccard", and a maximum distance between connected nodes of `0.3`.


```r
ig <- make_network(enterotype, max.dist = 0.3)
```


# Use the plot_network function
Because we want to use the enterotype designations as a plot feature in these plots, we need to remove the 9 samples for which no enterotype designation was assigned (this will save us the hassle of some pesky warning messages, but everything still works; the offending samples are anyway omitted).


```r
enterotype = subset_samples(enterotype, !is.na(Enterotype))
```


Now plot this network representation with the default settings.


```r
plot_network(ig, enterotype)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

The previous graphic displayed some interesting structure, with a major subgraph comprising a majority of samples. Furthermore, there seemed to be a correlation in the sample naming scheme and position within the network. Instead of trying to read all of the sample names to understand the pattern, let's map some of the sample variables onto this graphic as color and shape:


```r
plot_network(ig, enterotype, color = "SeqTech", shape = "Enterotype", line_weight = 0.4, 
    label = NULL)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


In the previous examples, the choice of maximum-distance and distance method were informed, but arbitrary. Let's see what happens when the maximum distance is lowered, decreasing the number of edges in the network


```r
ig <- make_network(enterotype, max.dist = 0.2)
plot_network(ig, enterotype, color = "SeqTech", shape = "Enterotype", line_weight = 0.4, 
    label = NULL)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Let's repeat the previous exercise, but replace the Jaccard (default) distance  method with Bray-Curtis


```r
ig <- make_network(enterotype, dist.fun = "bray", max.dist = 0.3)
plot_network(ig, enterotype, color = "SeqTech", shape = "Enterotype", line_weight = 0.4, 
    label = NULL)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

			


---

### Other tutorial pages for the phyloseq package:

#### [distance](distance.html)

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


