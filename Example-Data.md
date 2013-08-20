
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# Example Data for phyloseq


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.21'
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



# Example Data
## Included Example Data
There are multiple example data sets included in phyloseq. Many are from published investigations and include documentation with a summary and references, as well as some example code representing some aspect of analysis available in phyloseq.

To load example data into the working environment, use the `data()` command:


```r
data(GlobalPatterns)
data(esophagus)
data(enterotype)
data(soilrep)
```


## Example Data Documentation
In the package index, go to the names beginning with "data-" to see the documentation of currently available example datasets.

Try for example


```r
`?`(GlobalPatterns)
```


to access the documentation for the so-called "GlobalPatterns" dataset.

## Example Script Using Example Data
You can also try the examples included with the example data documentation (as well as examples for functions/methods) using the standard `example` command in R -- in this case the examples for the `enterotype` dataset.


```r
example(enterotype, ask = FALSE)
```

```
## 
## entrty> # Try simple network-analysis plot
## entrty> data(enterotype)
## 
## entrty> ig <- make_network(enterotype, "samples", max.dist=0.3)
## 
## entrty> plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
```

![plot of chunk run-examples](figure/run-examples1.png) 

```
## 
## entrty> #
## entrty> # Filter samples that don't have Enterotype
## entrty> x <- subset_samples(enterotype, !is.na(Enterotype))
## 
## entrty> #
## entrty> # Alternatively. . .
## entrty> ent.cca <- ordinate(x ~ Enterotype, "CCA")
## 
## entrty> plot_ordination(x, ent.cca, color="Enterotype")
```

![plot of chunk run-examples](figure/run-examples2.png) 

```
## 
## entrty> plot_ordination(x, ent.cca, "biplot")
```

![plot of chunk run-examples](figure/run-examples3.png) 

```
## 
## entrty> plot_ordination(x, ent.cca, "split", color="Enterotype")
```

![plot of chunk run-examples](figure/run-examples4.png) 

```
## 
## entrty> #
## entrty> # # multiple testing of genera correlating with enterotype 2
## entrty> # mt(x, data.frame(sample_data(x))[, "Enterotype"]==2)
## entrty> # # Should return a data.frame, with the following head()
## entrty> # # # # # index     teststat   rawp   adjp plower
## entrty> # # # Prevotella                      207 11.469961374 0.0001 0.0088 0.0001
## entrty> # # # Bacteroides                     203 -9.015717540 0.0001 0.0088 0.0001
## entrty> # # # Holdemania                      201 -5.810081084 0.0001 0.0088 0.0001
## entrty> # # # Acetivibrio                     156 -5.246137207 0.0001 0.0088 0.0001
## entrty> 
## entrty> 
## entrty>
```




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

