# phyloseq

See alternate documentation in the form of a [phyloseq wiki](https://github.com/joey711/phyloseq/wiki).

[A workshop session](https://secure.bioconductor.org/BioC2012/labs.php) regarding the use of [the phyloseq package](http://joey711.github.com/phyloseq/) in amplicon sequencing analysis will be provided at the

[Bioconductor Workshop (BioC 2012)](https://secure.bioconductor.org/BioC2012/)

on July 24-25, 2012 

At the Fred Hutchinson Cancer Research Center - Seattle, WA

## Description

The analysis of microbiological communities brings many challenges: the integration of many different types of data with methods from ecology, genetics, phylogenetics, network analysis, visualization and testing. The data itself may originate from widely different sources, such as the microbiomes of humans, soils, surface and ocean waters, wastewater treatment plants, industrial facilities, and so on; and as a result, these varied sample types may have very different forms and scales of related data that is extremely dependent upon the experiment and its question(s). The phyloseq package is a tool to import, store, analyze, and graphically display complex phylogenetic sequencing data that has already been clustered into Operational Taxonomic Units (OTUs), especially when there is associated sample data, phylogenetic tree, and/or taxonomic assignment of the OTUs. This package leverages many of the tools available in R for ecology and phylogenetic analysis (vegan, ade4, ape, picante), while also using advanced/flexible graphic systems (ggplot2) to easily produce publication-quality graphics of complex phylogenetic data. phyloseq uses a specialized system of S4 classes to store all related phylogenetic sequencing data as single experiment-level object, making it easier to share data and reproduce analyses. In general, phyloseq seeks to facilitate the use of R for efficient interactive and reproducible analysis of OTU-clustered high-throughput phylogenetic sequencing data.

More concretely, `phyloseq` provides:

 * Import abundance and related data from popular OTU clustering pipelines:
	- [QIIME](http://qiime.org/), [mothur](http://www.mothur.org/), [BIOM](http://www.qiime.org/svn_documentation/documentation/biom_format.html), [PyroTagger](http://pyrotagger.jgi-psf.org/cgi-bin/index.pl), [RDP](http://pyro.cme.msu.edu/). 
 	- Additional importers planned for [MG-RAST](http://metagenomics.anl.gov/), [CLoVR-16S](http://clovr.org/methods/clovr-16s/), [Genboree Microbiome](http://genboree.org/theCommons/projects/pub-gen-microbiome), [PANGEA](http://www.microgator.org/pangea/), and others.

 * Convenience analysis wrappers for common analysis tasks, e.g.
	- [distance](https://github.com/joey711/phyloseq/wiki/distance)  --> 44 distance methods supported (UniFrac, Jensen-Shannon, etc)
	- [ordinate](https://github.com/joey711/phyloseq/wiki/ordinate)  --> many supported methods, including constrained methods

 * Native, parallelized implementation of [UniFrac](https://github.com/joey711/phyloseq/wiki/Fast-Parallel-UniFrac) distance calculations.

 * Powerful, flexible custom plot methods using ggplot2 for rapid, convenient exploratory analysis.
	- [plot_heatmap](https://github.com/joey711/phyloseq/wiki/plot_heatmap)
	- [plot_tree](https://github.com/joey711/phyloseq/wiki/plot_tree)
	- [plot_ordination](https://github.com/joey711/phyloseq/wiki/plot_ordination)
	- [plot_network](https://github.com/joey711/phyloseq/wiki/plot_network)
	- [plot_richness_estimates](https://github.com/joey711/phyloseq/wiki/Graphics-Examples)
	- [plot_taxa_bar](https://github.com/joey711/phyloseq/wiki/plot_taxa_bar)

 * Multiple testing methods specific to high-throughput phylogenetic sequencing data.

 * Examples for analysis and graphics using real published data.

 * Some simple tree-based OTU clustering for building a combined OTU/tree data object "from scratch".

Additional features planned for the near-term:

 * Animated ordinations for time-series (and analogous) data.
	See beta-version support package: [animate.phyloseq](https://github.com/joey711/animate.phyloseq)

 * More advanced tree-based OTU clustering, leveraging tools already available in R.

 * Structures for additional data components, e.g. mass spec and expression data.

For news related to the latest development version, see:

[phyloseq NEWS](https://github.com/joey711/phyloseq/blob/master/inst/NEWS)

phyloseq is under active development, hosted and versioned at GitHub:

https://github.com/joey711/phyloseq

I make lots of effort to cite/attribute author contributions within official package documentation, citations, and anywhere else it is appropriate. Please feel free to fork and contribute!