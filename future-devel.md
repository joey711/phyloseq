
<link href="markdown.css" rel="stylesheet"></link>

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

1. **Big(ger)** Data
The firehose of new-gen sequencing data is making possible "big" datasets in this realm. For example, [the demo on importing the Human Microbiome Project data into R](http://joey711.github.com/phyloseq-demo/HMP_import_example.html) takes a considerable amount of time to run on a typical desktop/laptop and may push some less powerful machines to their limit. And that's just data import. We are considering some of the best approaches to help a tool like phyloseq address computational issues that are arising from dealing with this data of this size, without compromising some of the other features (interactivity, reproducibility, connection with existing R tools). Some promising tools already available for R that might help include:

##### Sparse matrix classes
Like in [the Matrix package](http://cran.r-project.org/web/packages/Matrix/index.html). This mainly applies to in-RAM computations, especially when the full matrix is actually needed. In principle, this might only apply be required to represent the preprocessed data, which could still be sparse.

##### Store full dataset in a database
Some candidate packages are:

- [the hdf5 package](http://cran.r-project.org/web/packages/hdf5/index.html)

- [the ff package](http://cran.r-project.org/web/packages/ff/index.html)



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

