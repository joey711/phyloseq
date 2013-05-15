
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>


Download microbio.me/qiime datasets
========================================================
[Originally hosted here](http://joey711.github.io/phyloseq/download-microbio.me.html)

Here is an example accessing microbiome datasets from a public repository, using entirely R code. I haven't yet figured out how to list the available studies at The [microbio.me/qiime](http://www.microbio.me/qiime/index.psp) from within R, but I have provided illustrate instructions for finding details about studies you might want to download and explore, and some example FTP addresses that I used after doing just that. Note that you likely need to create an account at The [microbio.me/qiime](http://www.microbio.me/qiime/index.psp) in order to explore the studies they have available for download.

## microbio.me/qiime --> *microbio_me_qiime* in phyloseq
 The `microbio_me_qiime` interface is a means to import data directly from the [microbio.me/qiime](http://www.microbio.me/qiime/index.psp) data repository with only a small quantity of effort.

There are datasets posted on [microbio.me/qiime](http://www.microbio.me/qiime/index.psp), and these usually have the "raw" OTU clustered data hosted at an FTP address. I'm not sure of a "nice" way to explor the details of the different studies using the FTP address directly, but I can post a few hand-picked datasets and import them all.

First, [create an account and login](http://www.microbio.me/qiime/index.psp). Here is what the frontpage should look like:

<img src="http://joey711.github.io/phyloseq/microbio-me-login.png" width="750px" />

Once you have logged-in, you will see a different screen. Click on "Get Study Summary and Raw Data", which will be a link near the top of the page (highlighted at the bottom of the following clipped page):

<img src="http://joey711.github.io/phyloseq/microbio-me-download-data.png" width="400px" />

And once you're at this "raw data" page, you will see a box below "Available Studies". The box is a large scrollable list, and the page below it will start off blank:

<img src="http://joey711.github.io/phyloseq/microbio-me-data-databox.png" width="250px" />

But once you click on something, details and additional links for that particular study will pop up. For example, I clicked on a version of the "Global Patterns" dataset in the box, which exposed a page below with further details about the study and links for data. If available, there is an FTP link called "Sequences, Mapping and OTU Table", which can be accessed directly from within R using the code provided farther down this page. The link itself is highlighted on this screenshot:

<img src="http://joey711.github.io/phyloseq/microbio-me-gpdata.png" width="250px" />


---

# Download, parse datasets

Load phyloseq package... and ggplot2 for plotting at the end...


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.15'
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
```


I previously defined a function here "on the fly" to download and parse data from the FTP server for `microbio.me`, and also download, unpack and import into R/[phyloseq](http://joey711.github.io/phyloseq/). A more complete version of this functionality is now included in phyloseq, and this page now demonstrates how that can work.

The main argument to this function identifies which study in the microbio.me/qiime database that you want to download and import. Here are several hand-picked examples using either just the study ID number.


```r
# Restroom Surfaces - Flores_restroom_surface_biogeography
restroom = microbio_me_qiime(1335)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# smokers - Charlson_cigarette_smokers
smokers = microbio_me_qiime(524)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# rs - resistant starches - Martinez_Resistant_starches
rs = microbio_me_qiime(495)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# abt - antibiotic timecourse - Relman_antibiotic_timeseries
abt = microbio_me_qiime(494)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# vagina - vaginal microbiome - Ravel_reproductive_women_vagina
vagina = microbio_me_qiime(509)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```

```r
# palms - Fierer_undergraduate_palms
palms = microbio_me_qiime(317)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```


Here is an example downloading a version of the Global Patterns dataset using the full URL (ftp) address:


```r
# # Global Patterns - CaporasoIlluminaPNAS2011_3prime
gpftp = "ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_721_split_library_seqs_and_mapping.tgz"
gp = microbio_me_qiime(gpftp)
```

```
## Found biom-format file, now parsing it... 
## Done parsing biom... 
## Importing Sample Metdadata from mapping file...
## Merging the imported objects... 
## Successfully merged, phyloseq-class created. 
##  Returning...
```


## Some minimal [preprocessing](http://joey711.github.io/phyloseq/preprocessing.html)
Remove "empty" OTUs and samples from each dataset. Note that I've chosen to keep only samples that have more than 100 total reads


```r
dlist = list(restroom, smokers, rs, abt, vagina, palms)
dlist = lapply(dlist, function(physeq) {
    physeq = prune_taxa(taxa_sums(physeq) > 0, physeq)
    physeq = prune_samples(sample_sums(physeq) > 100, physeq)
})
```


## Plots
And some plots, because humans like images and it also helps demonstrate that the data is now directly available for analysis.

### [Richness plots](http://joey711.github.io/phyloseq/plot_richness-examples.html)

```r
plot_richness(smokers, x = "AGE", color = "SEX", shape = "BODY_PRODUCT")
```

![plot of chunk richness-plots](figure/richness-plots1.png) 

```r
plot_richness(vagina, x = "NUGENT_SCORE", color = "ETHNICITY") + stat_smooth(method = lm)
```

![plot of chunk richness-plots](figure/richness-plots2.png) 

```r
plot_richness(vagina, x = "PH", color = "ETHNICITY") + stat_smooth(method = lm)
```

![plot of chunk richness-plots](figure/richness-plots3.png) 


### [Heatmap](http://joey711.github.io/phyloseq/plot_heatmap-examples.html)


```r
plot_heatmap(smokers)
```

![plot of chunk heatmap](figure/heatmap.png) 

