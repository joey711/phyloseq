
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>


# Gap Statistic
## How many clusters are there?

From [the clusGap documentation](http://stat.ethz.ch/R-manual/R-devel/library/cluster/html/clusGap.html): 
The `clusGap` function from [the cluster package](http://cran.r-project.org/web/packages/cluster/index.html) calculates a goodness of clustering measure, called [the “gap” statistic](www.stanford.edu/~hastie/Papers/gap.pdf). For each number of clusters `k`, it compares \log(W(k)) with E^*[\log(W(k))] where the latter is defined via bootstrapping, i.e. simulating from a reference distribution.

The following is an example performing the gap statistic on ordination results calculated using phyloseq tools, followed by an example of how a [ggplot](http://had.co.nz/ggplot2/)-based wrapper for this example might be included in [the phyloseq package](http://joey711.github.com/phyloseq/). 


### First perform an ordination

In this case, MDS on the Bray-Curtis distance.


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.18'
```

```r
library("cluster")
packageVersion("cluster")
```

```
## [1] '1.14.4'
```

```r
# Load data
data(enterotype)
# ordination
exord = ordinate(enterotype, method = "MDS", distance = "bray")
```



### Gap Statistic code


```r
pam1 = function(x, k) {
    list(cluster = pam(x, k, cluster.only = TRUE))
}
x = phyloseq:::scores(exord, display = "sites")
```

```
## Error: object 'scores' not found
```

```r
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN = pam1, K.max = 6, B = 50)
```

```
## Error: object 'x' not found
```

```r
gskmn
```

```
## Error: object 'gskmn' not found
```


Pretty straightforward. In case it is useful to see, this is what a wrapper-function might look like to "add-on" code to phyloseq.


```r
gap_statistic_ordination = function(ord, FUNcluster, type = "sites", K.max = 6, 
    axes = c(1:2), B = 500, verbose = interactive(), ...) {
    require("cluster")
    # If 'pam1' was chosen, use this internally defined call to pam
    if (FUNcluster == "pam1") {
        FUNcluster = function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
    }
    # Use the scores function to get the ordination coordinates
    x = phyloseq:::scores(ord, display = type)
    # If axes not explicitly defined (NULL), then use all of them
    if (is.null(axes)) {
        axes = 1:ncol(x)
    }
    # Finally, perform, and return, the gap statistic calculation using
    # cluster::clusGap
    clusGap(x[, axes], FUN = FUNcluster, K.max = K.max, B = B, verbose = verbose, 
        ...)
}
```


Define a plot method for results

```r
plot_clusgap = function(clusgap, title = "Gap Statistic calculation results") {
    require("ggplot2")
    gstab = data.frame(clusgap$Tab, k = 1:nrow(clusgap$Tab))
    p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
    p = p + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim))
    p = p + ggtitle(title)
    return(p)
}
```


Now try out this function. Should work on ordination classes recognized by `scores` function, and provide a [ggplot](http://had.co.nz/ggplot2/) graphic instead of a base graphic. (Special Note: the phyloseq-defined `scores` extensions are not exported as regular functions to avoid conflict, so phyloseq-defined `scores` extensions can only be accessed with the `phyloseq:::` namespace prefix in front.)


```r
gs = gap_statistic_ordination(exord, "pam1", B = 50, verbose = FALSE)
```

```
## Error: object 'scores' not found
```

```r
print(gs, method = "Tibs2001SEmax")
```

```
## Error: object 'gs' not found
```

```r
plot_clusgap(gs)
```

```
## Loading required package: ggplot2
```

```
## Error: object 'gs' not found
```


Base graphics plotting, for comparison.


```r
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
```

```
## Error: object 'gs' not found
```

```r
mtext("k = 2 is best ... but  k = 3  pretty close")
```

```
## Error: plot.new has not been called yet
```



---

## Other tutorial pages for the phyloseq package:

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



---

## Suggestions from Users/Developers

Don't be afraid to post feedback / needs on [the phyloseq issues tracker](https://github.com/joey711/phyloseq/issues):

https://github.com/joey711/phyloseq/issues
