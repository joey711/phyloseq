
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>


Functions for Accessing and (Pre)Processing Data
========================================================


```r
library(phyloseq)
```

For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line.


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.5.15"
```

Load the `GlobalPatterns` dataset, included with the phyloseq package.

```r
data(GlobalPatterns)
```


## Accessors
Components of a phyloseq object, like the OTU Table, can be accessed by special accessor functions, or ``accessors'', which return specific information about phylogenetic sequencing data, if present. These accessor functions are available for direct interaction by users and dependent functions/packages.


```r
GlobalPatterns
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

```r
ntaxa(GlobalPatterns)
```

```
## [1] 19216
```

```r
nsamples(GlobalPatterns)
```

```
## [1] 26
```

```r
sample_names(GlobalPatterns)[1:5]
```

```
## [1] "CL3"     "CC1"     "SV1"     "M31Fcsw" "M11Fcsw"
```

```r
rank_names(GlobalPatterns)
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

```r
sample_variables(GlobalPatterns)
```

```
## [1] "X.SampleID"               "Primer"                  
## [3] "Final_Barcode"            "Barcode_truncated_plus_T"
## [5] "Barcode_full_length"      "SampleType"              
## [7] "Description"
```

```r
otu_table(GlobalPatterns)[1:5, 1:5]
```

```
## OTU Table:          [5 taxa and 5 samples]
##                      taxa are rows
##        CL3 CC1 SV1 M31Fcsw M11Fcsw
## 549322   0   0   0       0       0
## 522457   0   0   0       0       0
## 951      0   0   0       0       0
## 244423   0   0   0       0       0
## 586076   0   0   0       0       0
```

```r
tax_table(GlobalPatterns)[1:5, 1:4]
```

```
## Taxonomy Table:     [5 taxa by 4 taxonomic ranks]:
##        Kingdom   Phylum          Class          Order         
## 549322 "Archaea" "Crenarchaeota" "Thermoprotei" NA            
## 522457 "Archaea" "Crenarchaeota" "Thermoprotei" NA            
## 951    "Archaea" "Crenarchaeota" "Thermoprotei" "Sulfolobales"
## 244423 "Archaea" "Crenarchaeota" "Sd-NA"        NA            
## 586076 "Archaea" "Crenarchaeota" "Sd-NA"        NA
```

```r
phy_tree(GlobalPatterns)
```

```
## 
## Phylogenetic tree with 19216 tips and 19215 internal nodes.
## 
## Tip labels:
## 	549322, 522457, 951, 244423, 586076, 246140, ...
## Node labels:
## 	, 0.858.4, 1.000.154, 0.764.3, 0.995.2, 1.000.2, ...
## 
## Rooted; includes branch lengths.
```

```r
taxa_names(GlobalPatterns)[1:10]
```

```
##  [1] "549322" "522457" "951"    "244423" "586076" "246140" "143239"
##  [8] "244960" "255340" "144887"
```

```r
myTaxa = taxa_names(GlobalPatterns)[1:10]
plot(phy_tree(prune_taxa(myTaxa, GlobalPatterns)))
```

![plot of chunk gp-sample-names](figure/gp-sample-names.png) 


## Preprocessing
The phyloseqBase package also includes functions for filtering, subsetting, and merging abundance data. Filtering in phyloseq is designed in a modular fashion similar to the approach in the genefilter package. This includes the `prune_taxa` and `prune_samples` methods for directly removing unwanted indices, as well as the `filterfun_sample` and `genefilter_sample` functions for building arbitrarily complex sample-wise filtering criteria, and the `filter_taxa` function for taxa-wise filtering. In the following example, the `GlobalPatterns` data is first transformed to relative abundance, creating the new `GPr` object, which is then filtered such that all OTUs with a variance greater than 10^-5 are removed. This results in a highly-subsetted object containing just 177 of the original ~19000 OTUs (`GPfr` below). Note that in both lines we have provided a custom function for transformation and filtering, respectively.


```r
GPr = transform_sample_counts(GlobalPatterns, function(x) x/sum(x))
GPfr = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)
```


The subsetting methods `prune_taxa` and `prune_samples` are for cases where the complete subset of desired OTUs or samples is directly available. Alternatively, the `subset_taxa` and `subset_samples` functions are for subsetting based on auxiliary data contained in the Taxonomy Table or Sample Data components, respectively. These functions are analogous to the `subset` function in core R, in which the initial data argument is followed by an arbitrary logical expression that indicates elements or rows to keep. Thus, entire experiment-level data objects can be subset according to conditional expressions regarding the auxiliary data. For example, the following code will first assign to `GP.chl` the subset of the GlobalPatterns dataset that are part of the Chlamydiae phylum, and then remove samples with less than 20 total reads.


```r
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
GP.chl = prune_samples(sampleSums(GP.chl) >= 20, GP.chl)
```


Merging methods include `merge_taxa` and `merge_samples`, intended for merging specific OTUs or samples, respectively.  There is also the `merge_phyloseq` function for a complete merge of two or more phyloseq-objects (or a phyloseq-object and one or more separate components). For example, the following code merges the first 5 OTUs in the Chlamydiae-only dataset.


```r
GP.chl.merged = merge_taxa(GP.chl, taxa_names(GP.chl)[1:5])
```


Building on the `merge_taxa` methods, the phyloseq-package includes the agglomeration functions, `tip_glom` and `tax_glom`, for merging all OTUs in an experiment that are similar beyond a phylogenetic or taxonomic threshold, respectively. The following code demonstrates how to agglomerate the "Bacteroidetes-only" dataset (called `gpsfb`) at the taxonomic rank of Family, and create an annotated tree of the result.


```r
gpsfbg = tax_glom(gpsfb, "Family")
plot_tree(gpsfbg, color = "SampleType", shape = "Class", size = "abundance")
```


For transforming abundance values by an arbitrary R function, phyloseqBase includes the `transform_sample_counts` function. It takes as arguments a phyloseq-object and an R function, and returns a phyloseq-object in which the abundance values have been transformed, sample-wise, according to the transformations specified by the function. For example, the following command transforms `GP.chl` abundance counts to fractional abundance.


```r
transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU))
```


Finally, the following is the remaining set of preprocessing steps that was applied to the GlobalPatterns OTU counts prior to creating the figures in the main phyloseq manuscript.

Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.

```r
GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2 * length(x)), 
    TRUE)
```


Define a human versus non-human categorical variable, and add this new variable to sample data:

```r
sample_data(GP)$human = factor(get_variable(GP, "SampleType") %in% c("Feces", 
    "Mock", "Skin", "Tongue"))
```


Standardize abundances to the median sequencing depth

```r
total = median(sample_sums(GP))
standf = function(x, t = total) round(t * (x/sum(x)))
gps = transform_sample_counts(GP, standf)
```


Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation

```r
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3, TRUE)
```


Subset the data to  Bacteroidetes, used in some plots

```r
gpsfb = subset_taxa(gpsf, Phylum == "Bacteroidetes")
```


## graphic summary
Now let's summarize this slice of the data with some graphics.


```r
title = "plot_bar; Bacteroidetes-only"
plot_bar(gpsfb, "SampleType", "Abundance", title = title)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) 

```r
plot_bar(gpsfb, "SampleType", "Abundance", "Family", title = title)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 



```r
plot_bar(gpsfb, "Family", "Abundance", "Family", title = title, facet_grid = "SampleType~.")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 




---

### Other tutorial pages for the phyloseq package:

#### [distance](distance.html)

#### [download-microbio.me-frag](download-microbio.me-frag.html)

#### [download-microbio.me](download-microbio.me.html)

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

#### [subset_ord_plot-examples](subset_ord_plot-examples.html)

#### [tutorials-index](tutorials-index.html)


