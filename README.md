# phyloseq

There are already several ecology and phylogenetic packages available in R, including the adephylo, vegan, ade4, picante, ape, phangorn, phylobase, and OTUbase packages. These can already take advantage of many of the powerful statistical and graphics tools available in R. However, at present a user must devise their own methods for parsing the output of their favorite OTU clustering application, and, as a consequence, there is also no standard within Bioconductor (or R generally) for storing or sharing the suite of related data objects that describe a phylogenetic sequencing project. The phyloseq package seeks to address these issues by providing a related set of S4 classes that internally manage the handling tasks associated with organizing, linking, storing, filtering and analyzing phylogenetic sequencing data. It further takes advantage of this structure by providing easy-to-use, but powerful, wrappers for analysis and publication-quality graphics. 'phyloseq' additionally provides some convenience wrappers for input from common clustering applications, common analysis pipelines, and native implementation of methods that are not available in other R packages.

More concretely, `phyloseq` provides:

 * Import abundance and related data from popular OTU clustering pipelines, like QIIME and mothur
	- BIOM, pyrotagger, RDP also supported.

 * Convenience analysis wrappers for common analysis tasks, including ordination.

 * Native, parallelized implementation of UniFrac distance calculations.

 * Some simple tree-based OTU clustering for building a combined OTU/tree data object "from scratch".

 * Powerful, flexible custom plot methods using ggplot2 for rapid, convenient exploratory analysis.
	- plot_sample_network
	- plot_richness_estimates
	- plot_ordination
	- plot_taxa_bar
	- plot_tree_phyloseq

 * Multiple testing methods specific to high-throughput phylogenetic sequencing data.

 * Examples for analysis and graphics using real published data.

Additional features planned for the near-term:

 * Animated ordinations for time-series (and analogous) data.
	See beta-version support package: https://github.com/joey711/animate.phyloseq

 * Structures for additional data components, e.g. mass spec and expression data.

For news related to the latest development version, see:

https://github.com/joey711/phyloseq/blob/master/inst/NEWS

phyloseq is under active development, which is hosted and versioned at github:

https://github.com/joey711/phyloseq

Feel free to fork and contribute!