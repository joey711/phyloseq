
<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>

Importing phyloseq Data
========================================================

The custom functions that read external data files and return an instance of the phyloseq-class are called "importers". Validity and coherency between data components are checked by the phyloseq-class constructor,
	phyloseq()
which is invoked internally by the importers, and is also the suggested function for creating a phyloseq object from "manually" imported data. The component indices representing OTUs or samples are checked for intersecting indices, and trimmed/reordered such that all available (non-`NULL`) component data describe exactly the same OTUs and samples, in the same order. 


```r
library(phyloseq)
```


For completeness, here is the version number of phyloseq used to build this instance of the tutorial -- and also how you can check your own current version from the command line. If your version is lower than the one shown here AND you are having trouble with the import examples, try updating your phyloseq version to the latest development version available from GitHub, using [the phyloseq installation tutorial](http://joey711.github.com/phyloseq/install).


```r
packageDescription("phyloseq")$Version
```

```
## [1] "1.5.4"
```


## Currently available import functions

See `?import` after phyloseq has been loaded (`library("phyloseq")`), to get an overview of available import functions and documentation links to their specific doc pages, or see below for examples using some of the more popular importers.

---
### import_biom

Newer versions of [QIIME](http://www.qiime.org/) produce a more-comprehensive and formally-defined JSON file format, called [biom file format](http://biom-format.org/):

"The biom file format (canonically pronounced ‘biome’) is designed to be a general-use format for representing counts of observations in one or more biological samples. BIOM is a recognized standard for the Earth Microbiome Project and is a Genomics Standards Consortium candidate project."

http://biom-format.org/

The phyloseq package includes small examples of biom files with different levels and organization of data. The following shows how to import each of the four main types of biom files (in practice, you don't need to know which type your file is, only that it is a biom file). In addition, the `import_biom` function allows you to simultaneously import an associated phylogenetic tree file and reference sequence file (e.g. fasta).

First, define the file paths. In this case, this will be within the phyloseq package, so we use special features of the `system.file` command to get the paths. This should also work on your system if you have phyloseq installed, regardless of your Operating System.


```r
rich_dense_biom = system.file("extdata", "rich_dense_otu_table.biom", package = "phyloseq")
rich_sparse_biom = system.file("extdata", "rich_sparse_otu_table.biom", package = "phyloseq")
min_dense_biom = system.file("extdata", "min_dense_otu_table.biom", package = "phyloseq")
min_sparse_biom = system.file("extdata", "min_sparse_otu_table.biom", package = "phyloseq")
treefilename = system.file("extdata", "biom-tree.phy", package = "phyloseq")
refseqfilename = system.file("extdata", "biom-refseq.fasta", package = "phyloseq")
```


Now that we've defined the file paths, let's use these as argument to the `import_biom` function. Note that the tree and reference sequence files are both suitable for any of the example biom files, which is why we only need one path for each. In practice, you will be specifying a path to a sequence or tree file that matches the rest of your data (include tree tip names and sequence headers)

```r
import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]
```

```r
import_biom(rich_sparse_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]
```

```r
import_biom(min_dense_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]
```

```r
import_biom(min_sparse_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]
```


Example code for importing large file with parallel backend

```r
library("doParallel")
registerDoParallel(cores = 6)
import_biom("my/file/path/file.biom", parseFunction = parse_taxonomy_greengenes, 
    parallel = TRUE)
```


In practice, you will store the result of your import as some variable name, like `myData`, and then use this data object in downstream data manipulations and analysis. For example,


```r
myData = import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)
plot_tree(myData, color = "Genus", shape = "BODY_SITE", size = "abundance")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) 

```r
plot_richness(myData, x = "BODY_SITE", color = "Description")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 

```r
plot_bar(myData, fill = "Genus")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-53.png) 

```r
refseq(myData)
```

```
##   A DNAStringSet instance of length 5
##     width seq                                          names               
## [1]   334 AACGTAGGTCACAAGCGTTGT...TTCCGTGCCGGAGTTAACAC GG_OTU_1
## [2]   465 TACGTAGGGAGCAAGCGTTAT...CCTTACCAGGGCTTGACATA GG_OTU_2
## [3]   249 TACGTAGGGGGCAAGCGTTAT...GGCTCGAAAGCGTGGGGAGC GG_OTU_3
## [4]   453 TACGTATGGTGCAAGCGTTAT...AAGCAACGCGAAGAACCTTA GG_OTU_4
## [5]   178 AACGTAGGGTGCAAGCGTTGT...GGAATGCGTAGATATCGGGA GG_OTU_5
```



---
### import_qiime

QIIME produces several files that can be analyzed in the phyloseq-package, including especially an OTU file that typically contains both OTU-abundance and taxonomic identity information. The map-file is also an important input to QIIME that stores sample covariates, converted naturally to the sample_data-class component data type in the phyloseq-package. QIIME may also produce a phylogenetic tree with a tip for each OTU, which can also be imported by this function.

See [qiime.org](http://www.qiime.org/) for details on using QIIME. While there are many complex dependencies, QIIME can be downloaded as a pre-installed linux virtual machine that runs “off the shelf”.

The different files useful for import to phyloseq are not collocated in a typical run of the QIIME pipeline. See the basics phyloseq vignette for an example of where to find the relevant files in the output directory.


```r
otufile = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package = "phyloseq")
mapfile = system.file("extdata", "master_map.txt", package = "phyloseq")
trefile = system.file("extdata", "GP_tree_rand_short.newick.gz", package = "phyloseq")
rs_file = system.file("extdata", "qiime500-refseq.fasta", package = "phyloseq")
qiimedata = import_qiime(otufile, mapfile, trefile, rs_file)
```

```
## Processing map file...
## Processing otu/tax file...
## 
## Reading and parsing file in chunks ... Could take some time. Please be patient...
## 
## Building OTU Table in chunks. Each chunk is one dot.
## .Building Taxonomy Table...
## Processing phylogenetic tree...
## Processing Reference Sequences...
```

```r
print(qiimedata)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 500 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 500 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 500 tips and 499 internal nodes ]
## refseq()      DNAStringSet:      [ 500 reference sequences ]
```


So it has Let's try some quick graphics built from our newly-imported dataset, `qiimedata`.


```r
plot_bar(qiimedata, x = "SampleType", fill = "Phylum")
```

![plot of chunk import-qiime-graphics](figure/import-qiime-graphics1.png) 

```r
plot_heatmap(qiimedata, sample.label = "SampleType", species.label = "Phylum")
```

![plot of chunk import-qiime-graphics](figure/import-qiime-graphics2.png) 



---
### import_mothur
The open-source, platform-independent, locally-installed software package, "mothur"", can also process barcoded amplicon sequences and perform OTU-clustering, among other things. It is extensively documented on a wiki at [the mothur wiki](http://www.mothur.org/wiki/Main_Page).


```r
mothlist = system.file("extdata", "esophagus.fn.list.gz", package = "phyloseq")
mothgroup = system.file("extdata", "esophagus.good.groups.gz", package = "phyloseq")
mothtree = system.file("extdata", "esophagus.tree.gz", package = "phyloseq")
show_mothur_list_cutoffs(mothlist)
```

```
##  [1] "unique" "0.00"   "0.01"   "0.02"   "0.03"   "0.04"   "0.05"  
##  [8] "0.06"   "0.07"   "0.08"   "0.09"   "0.10"
```

```r
cutoff = "0.10"
x = import_mothur(mothlist, mothgroup, mothtree, cutoff)
x
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 58 taxa and 3 samples ]
## phy_tree()    Phylogenetic Tree: [ 58 tips and 57 internal nodes ]
```

```r
plot_tree(x, color = "samples")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```r
SDF = data.frame(samples = sample_names(x), row.names = sample_names(x))
sample_data(x) = sample_data(SDF)
plot_richness(x)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 


The class and data in the object returned by `import_mothur` depends on the  arguments. If the first three arguments are provided, then a phyloseq object should be returned containing both a tree and its associated OTU table. If only a list and group file are provided, then an "otu_table" object is returned. Similarly, if only a list and tree file are provided, then only a tree is returned ("phylo" class).

Returns just a tree

```r
x1 = import_mothur(mothlist, mothur_tree_file = mothtree, cutoff = "0.10")
x2 = import_mothur(mothlist, mothur_tree_file = mothtree, cutoff = "0.08")
plot(x1)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

Returns just an OTU table

```r
OTU = import_mothur(mothlist, mothgroup, cutoff = "0.08")
dim(OTU)
```

```
## [1] 64  3
```

```r
head(OTU)
```

```
## OTU Table:          [6 taxa and 3 samples]
##                      taxa are rows
##          B  C   D
## 65_4_15 15 14   2
## 59_8_22 23  2   2
## 59_7_6  37 41  18
## 59_5_19 14 42  10
## 59_2_6  52 41 124
## 9_4_6    2  2   2
```


Returns a list where each (outer) element represents an OTU, and is a vector of the sequencing reads that are clustered with that OTU.

```r
otulist = import_mothur(mothlist, cutoff = "0.08")
length(otulist)
```

```
## [1] 64
```

```r
ntaxa(OTU)
```

```
## [1] 64
```


Returns an error without a cutoff

```r
import_mothur(mothlist)
```

```
## Error: non-character argument
```


The list file is required. Import will fail with an error if it is not provided.

```r
import_mothur()
```

```
## you must provide the mothur_list_file argument
```

```
## Error: argument "mothur_list_file" is missing, with no default
```



---
### import_pyrotagger

PyroTagger is created and maintained by the [Joint Genome Institute](http://pyrotagger.jgi-psf.org/)

The typical output form PyroTagger is a spreadsheet format ".xls", which poses additional import challenges. However, virtually all spreadsheet applications support the ".xls" format, and can further export this file in a tab-delimited format. It is recommended that you convert the xls-file without any modification (as tempting as it might be once you have loaded it) into a tab-delimited text file. Deselect any options to encapsulate fields in quotes, as extra quotes around each cell's contents might cause problems during file processing. These quotes will also inflate the file-size, so leave them out as much as possible, while also resisting any temptation to modify the xls-file “by hand”.

A highly-functional and free spreadsheet application can be obtained as part of [the cross-platform OpenOffice suite](http://www.openoffice.org/), and works for the above required conversion.

It is regrettable that this importer does not take the xls-file directly as input. However, because of the moving-target nature of spreadsheet file formats, there is limited support for direct import of these formats into R. Rather than add to the dependency requirements of emphphyloseq and the relative support of these xls-support packages, it seems more efficient to choose an arbitrary delimited text format, and focus on the data structure in the PyroTagger output. This will be easier to support in the long-run.

For example, the path to a pyrotagger tab-delimited file might be saved as `pyrotagger_tab_file`, and can be imported using:


```r
import_pyrotagger_tab(pyrotagger_tab_file)
```



---
## Loading included data

See [the wiki page on included example data in phyloseq](https://github.com/joey711/phyloseq/wiki/Example-Data) for more details.

The `data` command in the R language loads pre-imported datasets that are included in packages. For example, the "Global Patterns" dataset can be loaded into the R workspace with the following command.


```r
data(GlobalPatterns)
```




---

### Other tutorial pages for the phyloseq package:

#### [distance](distance.html)

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

#### [rebuild-all-from-Rmarkdown.R](rebuild-all-html-from-Rmarkdown.R)

#### [subset_ord_plot-examples](subset_ord_plot-examples.html)

#### [tutorials-index](tutorials-index.html)

