
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>


# Merge
[The phyloseq project](http://joey711.github.com/phyloseq/) includes support for two complete different categories of **merging**. 

- **Merging the OTUs or samples in a phyloseq object**, based upon a taxonomic or sample variable:
`merge_samples()`
`merge_taxa()`

Merging OTU or sample indices based on variables in the data can be a useful means of reducing noise or excess features in an analysis or graphic. Some examples might be to merge the samples in a dataset that are from the same environment, or orthogonally, to merge OTUs that are from the same taxonomic genera. 

- **Merging two or more data objects** that come from the same experiment, so that their data becomes part of the same phyloseq object:
`merge_phyloseq()`

Merging separate data objects is especially useful for manually-imported data objects, especially when one of the data objects already has more than one component and so is a phyloseq-class. While the first category of merging functions is useful for direct manipulations of the data for analytical purposes, `merge_phyloseq` is a convenience/support tool to help get your data into the right format.

# Examples

## merge_samples

`merge_samples` can be very useful if you'd like to see what happens to an analysis if you remove the indivual effects between replicates or between samples from a particular explanatory variable. With the `merge_samples` function, the abundance values of merged samples are summed, so make sure to do any [preprocessing](http://joey711.github.com/phyloseq/preprocess) to account for differences in sequencing effort before merging (or you will be doing a sequencing-effort-weighted average).

Let's first remove unobserved OTUs (sum 0 across all samples), and add a human-associated variable with which to organize data later in the plots.


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.21'
```


Load data, remove empty samples, add a new `sample_data` variable to the dataset.

```r
data(GlobalPatterns)
GP = GlobalPatterns
GP = prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
humantypes = c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% humantypes
```


Now on to the merging examples

```r
mergedGP = merge_samples(GP, "SampleType")
SD = merge_samples(sample_data(GP), "SampleType")
print(SD[, "SampleType"])
```

```
##                    SampleType
## Feces                       1
## Freshwater                  2
## Freshwater (creek)          3
## Mock                        4
## Ocean                       5
## Sediment (estuary)          6
## Skin                        7
## Soil                        8
## Tongue                      9
```

```r
print(mergedGP)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 18988 taxa and 9 samples ]
## sample_data() Sample Data:       [ 9 samples by 8 sample variables ]
## tax_table()   Taxonomy Table:    [ 18988 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 18988 tips and 18987 internal nodes ]
```

```r
sample_names(GP)
```

```
##  [1] "CL3"      "CC1"      "SV1"      "M31Fcsw"  "M11Fcsw"  "M31Plmr" 
##  [7] "M11Plmr"  "F21Plmr"  "M31Tong"  "M11Tong"  "LMEpi24M" "SLEpi20M"
## [13] "AQC1cm"   "AQC4cm"   "AQC7cm"   "NP2"      "NP3"      "NP5"     
## [19] "TRRsed1"  "TRRsed2"  "TRRsed3"  "TS28"     "TS29"     "Even1"   
## [25] "Even2"    "Even3"
```

```r
sample_names(mergedGP)
```

```
## [1] "Feces"              "Freshwater"         "Freshwater (creek)"
## [4] "Mock"               "Ocean"              "Sediment (estuary)"
## [7] "Skin"               "Soil"               "Tongue"
```

```r
identical(SD, sample_data(mergedGP))
```

```
## [1] TRUE
```


As emphasized earlier, the OTU abundances of merged samples are summed. Let's investigate this ourselves looking at just the top10 most abundance OTUs.


```r
OTUnames10 = names(sort(taxa_sums(GP), TRUE)[1:10])
GP10 = prune_taxa(OTUnames10, GP)
mGP10 = prune_taxa(OTUnames10, mergedGP)
ocean_samples = sample_names(subset(sample_data(GP), SampleType == "Ocean"))
print(ocean_samples)
```

```
## [1] "NP2" "NP3" "NP5"
```

```r
otu_table(GP10)[, ocean_samples]
```

```
## OTU Table:          [10 taxa and 3 samples]
##                      taxa are rows
##         NP2   NP3   NP5
## 329744   91   126   120
## 317182 3148 12370 63084
## 549656 5045 10713  1784
## 279599  113   114   126
## 360229   16    83   786
## 94166    49   128   709
## 550960   11    86    65
## 158660   13    39    28
## 331820   24   101   105
## 189047    4    33    29
```

```r
rowSums(otu_table(GP10)[, ocean_samples])
```

```
## 329744 317182 549656 279599 360229  94166 550960 158660 331820 189047 
##    337  78602  17542    353    885    886    162     80    230     66
```

```r
otu_table(mGP10)["Ocean", ]
```

```
## OTU Table:          [10 taxa and 1 samples]
##                      taxa are columns
##       329744 317182 549656 279599 360229 94166 550960 158660 331820 189047
## Ocean    337  78602  17542    353    885   886    162     80    230     66
```


Let's look at the merge graphically between two [richness estimate summary plots](http://joey711.github.com/phyloseq/plot_richness-examples). Note that we are ignoring our own advice and not-performing any preprocessing before this merge. Since these are richness estimates we are probably safe. 


```r
plot_richness(GP, "human", "SampleType", title = "unmerged")
```

```
## Warning: phyloseq::estimate_richness: Warning in fisher.alpha(). See
## `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


The merge can do some weird things to sample variables. Let's re-add these variables to the `sample_data` before we plot.


```r
sample_data(mergedGP)$SampleType = sample_names(mergedGP)
sample_data(mergedGP)$human = sample_names(mergedGP) %in% humantypes
plot_richness(mergedGP, "human", "SampleType", title = "merged")
```

```
## Warning: phyloseq::estimate_richness: Warning in fisher.alpha(). See
## `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Perhaps not surprisingly, when we combine the abundances of non-replicate samples from the same environment, the estimates of absolute richness increase for each environment. More to the example, however, the new merged plot is easier to read and interpret, which is one of the reasons one might use the `merge_samples` function.

## merge_taxa

One of the sources of "noise" in a microbial census dataset is a fine-scaled definition for OTU that blurs a pattern that might be otherwise evident were we to consider a higher taxonomic rank. The best way to deal with this is to use the agglomeration functions `tip_glom` or `tax_glom` that merge similar OTUs based on a phylogenetic or taxonomic threshold, respectively. However, for either function to work, it must be capable of merging two or more OTUs that have been deemed "equivalent". 

The following is an example of `merge_taxa` in action, using a tree graphic to display the before and after of the merge, and also show how the merge affects not just the OTU table, but all data components with OTU indices.

We will use a dataset that is included with the picante package, so a few steps must be taken to get the data into an optimal format for phyloseq functions.


```r
data(phylocom)
tree = phylocom$phylo
otu = otu_table(phylocom$sample, taxa_are_rows = FALSE)
x0 = phyloseq(otu, tree)
SDF = data.frame(sample = sample_names(x0), row.names = sample_names(x0))
sample_data(x0) = sample_data(SDF)
```


Now that the original unmerged dataset has been combined into a phyloseq object, `x0`, we are free to make an interesting tree graphic using the `plot_tree` function.

```r
plot_tree(x0, color = "sample", size = "abundance", sizebase = 2, label.tips = "taxa_names")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


Now let's use `merge_taxa` to merge the first 8 OTUS of `x0` into one new OTU. By choosing `2` as the optional third argument to `merge_taxa`, we are choosing to essentially pile the data of these 8 OTUs into the index for the second OTU. The default is to use the first OTU of the second argument.

```r
x1 = merge_taxa(x0, taxa_names(x0)[1:8], 2)
plot_tree(x1, color = "sample", size = "abundance", sizebase = 2, label.tips = "taxa_names")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


---
## merge_phyloseq

As said earlier, `merge_phyloseq` is a convenience/support tool to help get your data into the right format. Here is an example in which we extract components from an example dataset, and then build them back up to the original form using `merge_phyloseq` along the way.

Let's split apart the "Global Patterns" example dataset into some components


```r
data(GlobalPatterns)
tree = phy_tree(GlobalPatterns)
tax = tax_table(GlobalPatterns)
otu = otu_table(GlobalPatterns)
sam = sample_data(GlobalPatterns)
otutax = phyloseq(otu, tax)
otutax
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
```


As you can see, our new `otutax` object has just the OTU table and taxonomy table. Now let's use `merge_phyloseq` to build up the original `GlobalPatterns` object, and compare to make sure they are identical. Note how the arguments to `merge_phyloseq` are a mixture of multi-component (`otutax`) and single component objects.


```r
GP2 = merge_phyloseq(otutax, sam, tree)
identical(GP2, GlobalPatterns)
```

```
## [1] TRUE
```


### More complicated merges
The `merge_phyloseq` function will also work even with more than one multiple-component object.


```r
otusamtree = phyloseq(otu, sam, tree)
GP3 = merge_phyloseq(otusamtree, otutax)
GP3
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

```r
identical(GP3, GlobalPatterns)
```

```
## [1] FALSE
```


So this merge appears to have worked. The new object, `GP3`, looks similar to `GlobalPatterns` but is not identical. Why? Well, the assumption by `merge_phyloseq` is that you are atttempting to merge **separate** sources of abundance data, and so any portion of the OTU tables in the two phyloseq objects that have the same OTU indices are **summed** together, just like with `merge_taxa` earlier.

This example scenario was illustrative, but hopefully rare in practice. Nevertheless, just in case, an easy fix would be to extract the unique component of `otutax` and provide it to `merge_phyloseq`, instead of the entire phyloseq object. This amounts to a small intuitive modification to the previous `merge_phyloseq` command:


```r
GP4 = merge_phyloseq(otusamtree, tax_table(otutax))
GP4
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

```r
identical(GP4, GlobalPatterns)
```

```
## [1] TRUE
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

