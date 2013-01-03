
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

Importing phyloseq Data
========================================================

The custom functions that read external data files and return an instance of the phyloseq-class are called "importers". Validity and coherency between data components are checked by the phyloseq-class constructor,
	phyloseq()
which is invoked internally by the importers, and is also the suggested function for creating a phyloseq object from "manually" imported data. The component indices representing OTUs or samples are checked for intersecting indices, and trimmed/reordered such that all available (non-`NULL`) component data describe exactly the same OTUs and samples, in the same order. 


```r
library(phyloseq)
```



## Currently available import functions

See `?import` after phyloseq has been loaded (`library("phyloseq")`), to get an overview of available import functions, or see below for some of the more popular importers.


### import_biom


```r
# An included example of a rich dense biom file
rich_dense_biom = system.file("extdata", "rich_dense_otu_table.biom", package = "phyloseq")
import_biom(rich_dense_biom, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [5 taxa and 6 samples]
##                      taxa are rows
## Sample Data:         [6 samples by 4 sample variables]:
## Taxonomy Table:     [5 taxa by 7 taxonomic ranks]:
```

```r
# An included example of a sparse dense biom file
rich_sparse_biom = system.file("extdata", "rich_sparse_otu_table.biom", package = "phyloseq")
import_biom(rich_sparse_biom, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [5 taxa and 6 samples]
##                      taxa are rows
## Sample Data:         [6 samples by 4 sample variables]:
## Taxonomy Table:     [5 taxa by 7 taxonomic ranks]:
```


Example code for importing large file with parallel backend

```r
library("doParallel")
registerDoParallel(cores = 6)
import_biom("my/file/path/file.biom", parseFunction = parse_taxonomy_greengenes, 
    parallel = TRUE)
```



### import_qiime

QIIME produces several files that can be analyzed in the phyloseq-package, including especially an OTU file that typically contains both OTU-abundance and taxonomic identity information. The map-file is also an important input to QIIME that stores sample covariates, converted naturally to the sample_data-class component data type in the phyloseq-package. QIIME may also produce a phylogenetic tree with a tip for each OTU, which can also be imported by this function.

See [qiime.org](http://www.qiime.org/) for details on using QIIME. While there are many complex dependencies, QIIME can be downloaded as a pre-installed linux virtual machine that runs “off the shelf”.

The different files useful for import to phyloseq are not collocated in a typical run of the QIIME pipeline. See the basics phyloseq vignette for an example of where to find the relevant files in the output directory.


```r
otufile = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package = "phyloseq")
mapfile = system.file("extdata", "master_map.txt", package = "phyloseq")
trefile = system.file("extdata", "GP_tree_rand_short.newick.gz", package = "phyloseq")
import_qiime(otufile, mapfile, trefile, showProgress = FALSE)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [500 taxa and 26 samples]
##                      taxa are rows
## Sample Data:         [26 samples by 7 sample variables]:
## Taxonomy Table:     [500 taxa by 7 taxonomic ranks]:
## Phylogenetic Tree:  [500 tips and 499 internal nodes]
##                      rooted
```



### import_mothur
See [the mothur wiki](http://www.mothur.org/wiki/Main_Page) for further details about using mothur.


```r
mothlist = system.file("extdata", "esophagus.fn.list.gz", package = "phyloseq")
mothgroup = system.file("extdata", "esophagus.good.groups.gz", package = "phyloseq")
mothtree = system.file("extdata", "esophagus.tree.gz", package = "phyloseq")
cutoff = "0.10"
import_mothur(mothlist, mothgroup, mothtree, cutoff)
```

```
## phyloseq-class experiment-level object
## OTU Table:          [58 taxa and 3 samples]
##                      taxa are rows
## Phylogenetic Tree:  [58 tips and 57 internal nodes]
##                      rooted
```


Will fail with an error if no list file provided.

```r
import_mothur()
```

```
## you must provide the mothur_list_file argument
```

```
## Error: argument "mothur_list_file" is missing, with no default
```



### import_pyrotagger

PyroTagger is created and maintained by the [Joint Genome Institute](http://pyrotagger.jgi-psf.org/)

The typical output form PyroTagger is a spreadsheet format ".xls", which poses additional import challenges. However, virtually all spreadsheet applications support the ".xls" format, and can further export this file in a tab-delimited format. It is recommended that you convert the xls-file without any modification (as tempting as it might be once you have loaded it) into a tab-delimited text file. Deselect any options to encapsulate fields in quotes, as extra quotes around each cell's contents might cause problems during file processing. These quotes will also inflate the file-size, so leave them out as much as possible, while also resisting any temptation to modify the xls-file “by hand”.

A highly-functional and free spreadsheet application can be obtained as part of [the cross-platform OpenOffice suite](http://www.openoffice.org/), and works for the above required conversion.

It is regrettable that this importer does not take the xls-file directly as input. However, because of the moving-target nature of spreadsheet file formats, there is limited support for direct import of these formats into R. Rather than add to the dependency requirements of emphphyloseq and the relative support of these xls-support packages, it seems more efficient to choose an arbitrary delimited text format, and focus on the data structure in the PyroTagger output. This will be easier to support in the long-run.

For example, the path to a pyrotagger tab-delimited file might be saved as `pyrotagger_tab_file`, and can be imported using:


```r
import_pyrotagger_tab(pyrotagger_tab_file)
```



## Loading included data

See [the wiki page on included example data in phyloseq](https://github.com/joey711/phyloseq/wiki/Example-Data) for more details.

The `data` command in the R language loads pre-imported datasets that are included in packages. For example, the "Global Patterns" dataset can be loaded into the R workspace with the following command.


```r
data(GlobalPatterns)
```


