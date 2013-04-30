
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# The distance function in phyloseq

The `distance` function takes a phyloseq-class object and method option, and returns a `dist`-class distance object suitable for certain ordination methods and other distance-based analyses. There are currently 44 explicitly supported method options in the phyloseq package, as well as user-provided arbitrary methods via an interface to `vegan::designdist`. For the complete list of currently supported options/arguments to the method parameter, type `distance("list")` in the command-line of your R session. Only sample-wise distances are currently supported (the type argument), but eventually OTU-wise (e.g. species) distances will be supported as well.

See the in-package documentation of `distance` for further details:


```r
`?`(distance)
```


## Usage


```r
distance(physeq, method = "unifrac", type = "samples", ...)
```


## Example: "Enterotypes" dataset using many different methods
Because the `distance()` function organizes distance calculations into one function, it is relatively straightforward to calculate all supported distance methods and investigate the results. The following code will perform such a loop on the "Enterotypes" dataset, perform multi-dimensional scaling (a.k.a. principle coordinates analysis), and plot the first two axes, shading and shaping the points in each plot according to sequencing technology and assigned "Enterotype" label.

Note that we have omitted the options that require a phylogenetic tree because the `"enterotype"` example dataset currently included in the phyloseq-package does not have one.

Note that this may take a little while to run, depending on the size of your data set, but you may not be interested in all supported distances...


---

# Load phyloseq
Of course we need to start this tutorial by loading [the phyloseq package](http://joey711.github.com/phyloseq/). This assumes you have already [installed phyloseq](http://joey711.github.io/phyloseq/install).


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.8'
```

```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```

```r
library("plyr")
packageVersion("plyr")
```

```
## [1] '1.8'
```


ggplot2 package theme set. See [the ggplot2 online documentation](http://docs.ggplot2.org/current/) for further help.

```r
theme_set(theme_bw())
```



Load the enterotype data

```r
data(enterotype)
```


Some preliminary filtering. More advanced [preprocessing](http://joey711.github.com/phyloseq/preprocess) is recommended.

Remove the OTUs that included all unassigned sequences (`"-1"`)

```r
enterotype <- subset_species(enterotype, Genus != "-1")
```


The available distance methods coded in `distance`

```r
dist_methods <- unlist(distance("list"))
print(dist_methods)
```

```
##      UniFrac        DPCoA          JSD     vegdist1     vegdist2 
##    "unifrac"      "dpcoa"        "jsd"  "manhattan"  "euclidean" 
##     vegdist3     vegdist4     vegdist5     vegdist6     vegdist7 
##   "canberra"       "bray" "kulczynski"    "jaccard"      "gower" 
##     vegdist8     vegdist9    vegdist10    vegdist11    vegdist12 
##   "altGower"   "morisita"       "horn"  "mountford"       "raup" 
##    vegdist13    vegdist14    vegdist15   betadiver1   betadiver2 
##   "binomial"       "chao"        "cao"          "w"         "-1" 
##   betadiver3   betadiver4   betadiver5   betadiver6   betadiver7 
##          "c"         "wb"          "r"          "I"          "e" 
##   betadiver8   betadiver9  betadiver10  betadiver11  betadiver12 
##          "t"         "me"          "j"        "sor"          "m" 
##  betadiver13  betadiver14  betadiver15  betadiver16  betadiver17 
##         "-2"         "co"         "cc"          "g"         "-3" 
##  betadiver18  betadiver19  betadiver20  betadiver21  betadiver22 
##          "l"         "19"         "hk"        "rlb"        "sim" 
##  betadiver23  betadiver24        dist1        dist2        dist3 
##         "gl"          "z"    "maximum"     "binary"  "minkowski" 
##   designdist 
##        "ANY"
```


Remove the two distance-methods that require a tree, and the generic custom method that requires user-defined distance arguments.

```r
# These require tree
dist_methods[(1:2)]
```

```
##   UniFrac     DPCoA 
## "unifrac"   "dpcoa"
```

```r
# Remove them from the vector
dist_methods <- dist_methods[-(1:2)]
# This is the user-defined method:
dist_methods["designdist"]
```

```
## designdist 
##      "ANY"
```

```r
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods == "ANY")]
```


Loop through each distance method, save each plot to a list, called `plist`.

```r
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for (i in dist_methods) {
    # Calculate distance matrix
    iDist <- distance(enterotype, method = i)
    # Calculate ordination
    iMDS <- ordinate(enterotype, "MDS", distance = iDist)
    ## Make plot Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    # Create plot, store as temp variable, p
    p <- plot_ordination(enterotype, iMDS, color = "SeqTech", shape = "Enterotype")
    # Add title to each plot
    p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
    # Save the graphic to file.
    plist[[i]] = p
}
```



## Combine all results

Shade according to **sequencing technology**

```r
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = SeqTech, shape = Enterotype))
p = p + geom_point(size = 3, alpha = 0.5)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p
```

![plot of chunk distance-summary-plot](figure/distance-summary-plot.png) 


Shade according to assigned **enterotype**

```r
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = Enterotype, shape = SeqTech))
p = p + geom_point(size = 3, alpha = 0.5)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p
```

![plot of chunk distance-summary-plot-enterotype](figure/distance-summary-plot-enterotype.png) 



## Selected Results

The following are some selected examples among the created plots.

Jensen-Shannon Divergence

```r
print(plist[["jsd"]])
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


Jaccard

```r
print(plist[["jaccard"]])
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


Bray-Curtis

```r
print(plist[["bray"]])
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


Gower

```r
print(plist[["gower"]])
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


w

```r
print(plist[["w"]])
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 




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

#### [subset_ord_plot-examples](subset_ord_plot-examples.html)

#### [tutorials-index](tutorials-index.html)


