
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>


Subset points in an ordination plot
========================================================

The `subset_ord_plot` function is a "convenience function" intended to make it easier to retrieve a plot-derived `data.frame` with a subset of points according to a `threshold` and `method`. The meaning of the `threshold` depends upon the `method`.

Load the necessary packages and data.

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.7'
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



Some subsetting and light massaging of the `GlobalPatterns` dataset.

Clean zeros.

```r
GP <- GlobalPatterns
GP <- prune_species(taxa_sums(GP) > 0, GP)
```


Add `"human"` variable to `GP` to indicate human-associated samples.

```r
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", 
    "Skin", "Tongue")
```


Subset to just Bacteroidetes phylum. 

```r
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

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Re-draw this as topo without points, and facet

```r
p0 = ggplot(p1$data, p1$mapping) + geom_density2d() + facet_wrap(~Class)
p0
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Re-draw this but include points

```r
p1 = p1 + geom_density2d() + facet_wrap(~Class)
p1
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


Add a layer of a subset of species-points that are furthest from origin.

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "square"), size = 1)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-91.png) 

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "farthest"), size = 1)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-92.png) 

```r
p0 + geom_point(data = subset_ord_plot(p1, 0.7, "radial"), size = 1)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-93.png) 


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


