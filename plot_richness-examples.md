
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

plot_richness examples
========================================================

## Graphic Summary of Richness Estimates

As usual, we must start by loading the phyloseq package, and then the dataset, in this case `"GlobalPatterns"`.


```r
library("phyloseq")
data("GlobalPatterns")
```


For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.3.12"
```


Since we are interested in richness estimates, it is probably not a bad idea to prune OTUs that are not present in any of the samples (for some reason there are a few in `"GlobalPatterns"`) -- **BUT DON'T TRIM MORE THAN THAT!** I know it is tempting to trim noise right away, but many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.


```r
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
```


Here is the default graphic produced by the `plot_richness` function on the `GP` example dataset:


```r
plot_richness(GP)
```

![plot of chunk default](figure/default.png) 


We can specify a sample variable on which to group/organize samples along the horizontal (`x`) axis. An experimentally meaningful categorical variable is usually a good choice -- in this case, the `"SampleType"` variable works much better than attempting to interpret the sample names directly (as in the previous plot):


```r
plot_richness(GP, "SampleType")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
plot_richness(GP, x = "SampleType", color = "SampleType")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

We can also add the Shannon and Simpson alpha diversity indices by setting the `shsi` argument to `TRUE`. Naturally, the scales of these indices are very different from estimates for numbers of OTUs, and so the axes limits are allowed to differ in each panel when this option is `TRUE`.

```r
plot_richness(GP, x = "SampleType", color = "SampleType", shsi = TRUE)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


Now suppose we wanted to use an external variable in the plot that isn't in the `GP` dataset already -- for example, a logical that indicated whether or not the samples are human-associated. First, define this new variable, `human`, as a factor (other vectors could also work).


```r
# Define a human-associated versus non-human categorical variable:
sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", 
    "Skin", "Tongue")
```


Now tell `plot_richness` to map the new `human` variable on the horizontal axis, and shade the points in different color groups, according to which `"SampleType"` they belong.


```r
plot_richness(GP, "human", "SampleType")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

		
			

## Other key graphics functions in the phyloseq package:

### [plot_ordination](http://joey711.github.com/phyloseq/plot_ordination-examples)

### [plot_heatmap](http://joey711.github.com/phyloseq/plot_heatmap-examples)

### [plot_network](http://joey711.github.com/phyloseq/plot_network-examples)

### [plot_tree](http://joey711.github.com/phyloseq/plot_tree-examples)

### [plot_bar](http://joey711.github.com/phyloseq/plot_bar-examples)

### [plot_richness](http://joey711.github.com/phyloseq/plot_richness-examples)


