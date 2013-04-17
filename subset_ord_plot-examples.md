
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>


Subset points in an ordination plot
========================================================

The `subset_ord_plot` function is a "convenience function" intended to make it easier to retrieve a plot-derived `data.frame` with a subset of points according to a `threshold` and `method`. The meaning of the `threshold` depends upon the `method`.

Load the phyloseq package, and the GlobalPatterns example dataset


```r
library("phyloseq")
data("GlobalPatterns")
```

For completeness, here is the version number of phyloseq used to build this instance of the tutorial – and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.5.3"
```


Some subsetting and light massaging of the `GlobalPatterns` dataset.

```r
# Need to clean the zeros from GlobalPatterns:
GP <- GlobalPatterns
GP <- prune_species(taxa_sums(GP) > 0, GP)
# Add 'human' variable to GP
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", 
    "Skin", "Tongue")
# Subset to just Bacteroidetes
GP <- subset_taxa(GP, Phylum == "Bacteroidetes")
```


Perform a correspondence analysis and then create some plots to demonstrate using `subset_ord_plot`. We want to make species topo with a subset of points layered. Start by performing the correspondence analysis.


```r
gpca <- ordinate(GP, "CCA")
```


Now make a basic plot of just the species points in the correspondence analysis.

```r
p1 = plot_ordination(GP, gpca, "species", color = "Class")
p1
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Re-draw this as topo without points, and facet

```r
library("ggplot2")
p0 = ggplot(p1$data, p1$mapping) + geom_density2d() + facet_wrap(~Class)
p0
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Re-draw this but include points

```r
p1 = p1 + geom_density2d() + facet_wrap(~Class)
p1
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Add a layer of a subset of species-points that are furthest from origin.

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "square"), size = 1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) 

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "farthest"), size = 1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "radial"), size = 1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-83.png) 


Here is what the data retreived by `subset_ord_plot` actually looks like

```r
head(subset_ord_plot(p1, 0.7, "radial"))
```

```
##           CA1    CA2  Kingdom        Phylum           Class
## 554668 1.5118 1.0088 Bacteria Bacteroidetes Sphingobacteria
## 154451 1.5386 1.0803 Bacteria Bacteroidetes Sphingobacteria
## 244965 0.9385 0.3647 Bacteria Bacteroidetes Sphingobacteria
## 62487  1.5268 1.0448 Bacteria Bacteroidetes Sphingobacteria
## 139513 1.5244 1.0375 Bacteria Bacteroidetes Sphingobacteria
## 330672 0.8197 0.1491 Bacteria Bacteroidetes Sphingobacteria
##                     Order       Family    Genus             Species
## 554668 Sphingobacteriales Balneolaceae Balneola                <NA>
## 154451 Sphingobacteriales Balneolaceae Balneola Balneolaalkaliphila
## 244965 Sphingobacteriales Balneolaceae Balneola Balneolaalkaliphila
## 62487  Sphingobacteriales Balneolaceae Balneola                <NA>
## 139513 Sphingobacteriales Balneolaceae Balneola                <NA>
## 330672 Sphingobacteriales Balneolaceae     <NA>                <NA>
```


