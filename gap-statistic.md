
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
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```

```r
theme_set(theme_bw())
# Load data
data(enterotype)
# ordination
exord = ordinate(enterotype, method = "MDS", distance = "jsd")
```



### Gap Statistic code


```r
pam1 = function(x, k) {
    list(cluster = pam(x, k, cluster.only = TRUE))
}
x = phyloseq:::scores.pcoa(exord, display = "sites")
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN = pam1, K.max = 6, B = 50)
gskmn
```

```
## Clustering Gap statistic ["clusGap"].
## B=50 simulated reference sets, k = 1..6
##  --> Number of clusters (method 'firstSEmax', SE.factor=1): 4
##       logW E.logW    gap  SE.sim
## [1,] 2.996  3.117 0.1210 0.02110
## [2,] 2.210  2.771 0.5609 0.02413
## [3,] 1.922  2.587 0.6649 0.02467
## [4,] 1.686  2.412 0.7266 0.02602
## [5,] 1.601  2.283 0.6822 0.02622
## [6,] 1.481  2.178 0.6972 0.02884
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
    x = phyloseq:::scores.pcoa(ord, display = type)
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
print(gs, method = "Tibs2001SEmax")
```

```
## Clustering Gap statistic ["clusGap"].
## B=50 simulated reference sets, k = 1..6
##  --> Number of clusters (method 'Tibs2001SEmax', SE.factor=1): 4
##       logW E.logW    gap  SE.sim
## [1,] 2.996  3.111 0.1151 0.02088
## [2,] 2.210  2.770 0.5597 0.02211
## [3,] 1.922  2.581 0.6587 0.02575
## [4,] 1.686  2.415 0.7287 0.03023
## [5,] 1.601  2.284 0.6834 0.02097
## [6,] 1.481  2.181 0.7002 0.02161
```

```r
plot_clusgap(gs)
```

![plot of chunk gapstat-inphyloseq-example](figure/gapstat-inphyloseq-example.png) 


Base graphics plotting, for comparison.


```r
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("Looks like 4 clusters is best, with 3 and 5 close runners up.")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 



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
