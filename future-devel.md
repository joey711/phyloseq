
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

Future Development of phyloseq
========================================================

Additional development of phyloseq is ongoing, the details of which are documented heavily at [the phyloseq issues tracker](https://github.com/joey711/phyloseq/issues).

In broader strokes, the near-term plans include:

1. The "compartmentalization" of the data infrastructure portion of the phyloseq package into a separate Bioconductor package called [phyloseqBase](https://github.com/joey711/phyloseq/issues/102). This may reduce the number of dependencies implied when reusing the data tools of phyloseq for other Bioconductor/R packages.

1. The addition to phyloseqBase of interfaces to public repositories for integrated microbiome census data; for instance, [MG-RAST](http://metagenomics.anl.gov/) and [microbe.me](http://microbio.me/qiime).

1. More advanced tree-based OTU clustering, leveraging tools already available in R.

1. The development and calibration of additional tools in phyloseqBase for formal preprocessing, filtering, normalization, shrinkage and variance stabilization ([e.g. Allison 2006](http://www.nature.com/nrg/journal/v7/n1/full/nrg1749.html)) 

1. Additional structured wrapping tools to the [ade4](http://cran.r-project.org/web/packages/ade4/index.html) and [vegan](http://cran.r-project.org/web/packages/vegan/index.html) packages, with the most commonly-used tools given specific wrapping functions in phyloseq or phyloseqBase.

1. Structures for additional data components, e.g. mass spec and expression data.

1. Animated ordinations for time-series (and analogous) data. See beta-version support package: [animate.phyloseq](https://github.com/joey711/animate.phyloseq)

1. The compilation of a [Bioconductor data package](http://www.bioconductor.org/packages/release/data/experiment/) ("phyloseqData") that includes many key published datasets already imported as separate "phyloseq" instances and available through R's `data` interface.
