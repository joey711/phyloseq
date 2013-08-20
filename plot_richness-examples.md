
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

Alpha diversity graphics
========================================================

## plot_richness() examples
Using the `plot_richness` function.

Although the function name includes the word `richness`, which usually refers to the total number of species/OTUs in a sample or environment -- either observed or estimated -- this is actually a wrapper for all descriptions of [alpha diversity](http://en.wikipedia.org/wiki/Alpha_diversity). The name of this function may be changed in future versions to reflect this and avoid confusion.

## Graphic Summary of Alpha Diversity estimators

As usual, we must start by loading the phyloseq package, and then the dataset, in this case `"GlobalPatterns"`.


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.21'
```

```r
data("GlobalPatterns")
```


Some ggplot2 theming. First load [the ggplot2 package](http://docs.ggplot2.org/current/).


```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```

```r
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
    scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
    scale_fill_brewer(palette = palname, ...)
}
```



Since we are interested in alpha diversity, it is probably not a bad idea to prune OTUs that are not present in any of the samples (for some reason there are a few in `"GlobalPatterns"`) -- **BUT DON'T TRIM MORE THAN THAT!** I know it is tempting to trim noise right away, but many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.


```r
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
```


Here is the default graphic produced by the `plot_richness` function on the `GP` example dataset:


```r
plot_richness(GP)
```

```
## Warning: phyloseq::estimate_richness: Warning in fisher.alpha(). See
## `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution
```

![plot of chunk default](figure/default.png) 


Note that in this case, the Fisher calculation results in a warning (but still plots). We can avoid this by specifying a `measures` argument to `plot_richness`, which will include just the alpha-diversity measures that we want.


```r
plot_richness(GP, measures = c("Chao1", "Shannon"))
```

![plot of chunk two-measures-sample](figure/two-measures-sample.png) 


We can specify a sample variable on which to group/organize samples along the horizontal (`x`) axis. An experimentally meaningful categorical variable is usually a good choice -- in this case, the `"SampleType"` variable works much better than attempting to interpret the sample names directly (as in the previous plot):


```r
plot_richness(GP, x = "SampleType", measures = c("Chao1", "Shannon"))
```

![plot of chunk separate-sample](figure/separate-sample.png) 


Now suppose we wanted to use an external variable in the plot that isn't in the `GP` dataset already -- for example, a logical that indicated whether or not the samples are human-associated. First, define this new variable, `human`, as a factor (other vectors could also work; or other data you might have describing the samples).


```r
sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", 
    "Skin", "Tongue")
```


Now tell `plot_richness` to map the new `human` variable on the horizontal axis, and shade the points in different color groups, according to which `"SampleType"` they belong.


```r
plot_richness(GP, x = "human", color = "SampleType", measures = c("Chao1", "Shannon"))
```

![plot of chunk plot-human-1](figure/plot-human-1.png) 


We can merge samples that are from the environment (`SampleType`), and make the points bigger with a ggplot2 layer. First, merge the samples.


```r
GPst = merge_samples(GP, "SampleType")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(GPst)$SampleType <- factor(sample_names(GPst))
sample_data(GPst)$human <- as.logical(sample_data(GPst)$human)
```


Now we can plot this environment-merged version of the data. First store the default ggplot graphic as `p`, then add an additional `geom_point` layer with a large size and slight transparency. 


```r
p = plot_richness(GPst, x = "human", color = "SampleType", measures = c("Chao1", 
    "Shannon"))
p + geom_point(size = 5, alpha = 0.7)
```

![plot of chunk plot-human-2](figure/plot-human-2.png) 



### More details about ggplot2

For those interested in why this works so concisely (`p + geom_point(size=4, alpha=0.7)`), it is because the rest of the aesthetic mapping and data are contained in the ggplot object, `p`, and so is inherited in the call to the ggplot2 geometric object layer function, `geom_point`, by default since we didn't specify alternative `aes` or `data` arguments. Although we could have if we wanted to. This perhaps sounds more confusing than it is, and I find it easier to understand by inspecting the examples I've shown here. 

You'll also notice that the original smaller points are still on the plot. This is because they were the first layer, and our larger points are semi-transparent. I find this kind of distracting, and doesn't add any information or clarity. The good news is that layers can be removed from a ggplot object with standard list notation (using the dollar sign `$`).

First, check which lists are present in `p`.

```r
p$layers
```

```
## [[1]]
## geom_point: na.rm = TRUE 
## stat_identity:  
## position_identity: (width = NULL, height = NULL)
## 
## [[2]]
## mapping: ymax = value + se, ymin = value - se 
## geom_errorbar: width = 0.1 
## stat_identity:  
## position_identity: (width = NULL, height = NULL)
```


We can see that the first layer is the one specifying the original points, which are small. We can use negative indexing to "pop" it out, then add a new `geom_point` layer with larger point size (the following two lines).


```r
p$layers <- p$layers[-1]
p + geom_point(size = 5, alpha = 0.7)
```

![plot of chunk remove-layers](figure/remove-layers.png) 


			

---

### Other tutorial pages for the phyloseq package:

#### [distance](distance.html)

#### [download-microbio.me-frag](download-microbio.me-frag.html)

#### [download-microbio.me](download-microbio.me.html)

#### [Example-Data](Example-Data.html)

#### [future-devel](future-devel.html)

#### [gap-statistic](gap-statistic.html)

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

#### [subset_ord_plot-examples](subset_ord_plot-examples.html)

#### [tutorials-index](tutorials-index.html)


