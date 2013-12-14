
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# Importing phyloseq Data

The custom functions that read external data files and return an instance of the phyloseq-class are called *importers*. Validity and coherency between data components are checked by the phyloseq-class constructor, `phyloseq()` which is invoked internally by the importers, and is also the suggested function for creating a phyloseq object from [manually imported data](#manual). The component indices representing OTUs or samples are checked for intersecting indices, and trimmed/reordered such that all available (non-) component data describe exactly the same OTUs and samples, in the same order. 

See `?import` after phyloseq has been loaded (`library("phyloseq")`), to get an overview of available import functions and documentation links to their specific doc pages, or see below for examples using some of the more popular [import functions](#import_functions).

## [Create phyloseq Data Manually](#manual)

### [MG-RAST](#mgrast)

## <a name="import_functions"></a> Currently available import functions

### [microbio_me_qiime](#microbio_me)

### [import_biom](#import_biom)

### [import_qiime](#import_qiime)

### [import_mothur](#import_mothur)

### [import_pyrotagger](#import_pyrotagger)


---


---


---

---
## Begin examples, load requisite packages


```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.7.10'
```

```r
library("ggplot2")
packageVersion("ggplot2")
```

```
## [1] '0.9.3.1'
```


Define a default theme for ggplot graphics.

```r
theme_set(theme_bw())
```



---

### <a name="manual"></a> Create phyloseq Data Manually

There are lots of ways to get data related to a microbiome project, and not all of these will come from a popular server or workflow that is already supported in phyloseq. We also want to encourage users to create and share their own import code for special data formats. For these reasons especially, phyloseq provides tools for *constructing* phyloseq component data, and the experiment-level multi-component data object, the *phyloseq-class*. These are the same functions used internally by [the currently available importers](#import_functions).


```r
`?`(phyloseq)
`?`(otu_table)
`?`(sample_data)
`?`(tax_table)
```


<img src="http://www.plosone.org/article/info:doi/10.1371/journal.pone.0061217.g003/largerimage" width="750px" />

**If you can get the data into R, then you can get it "into" phyloseq.**

Constructors:

- `otu_table` - Works on any numeric `matrix`. You must also specify if the species are rows or columns
- `sample_data` - Works on any `data.frame`. The rownames must match the sample names in the `otu_table` if you plan to combine them as a phyloseq-object
- `tax_table` - Works on any character `matrix`. The rownames must match the OTU names (`taxa_names`) of the `otu_table` if you plan to combine it with a phyloseq-object.
- `phyloseq` - Takes as argument an `otu_table` and any unordered list of valid phyloseq components: `sample_data`, `tax_table`, `phylo`, or `XStringSet`. The tip labels of a phylo-object (tree) must match the OTU names of the `otu_table`, and similarly, the sequence names of an `XStringSet` object must match the OTU names of the `otu_table`. 
- `merge_phyloseq` - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.

**Note:** OTUs and samples are included in the combined object only if they are present in all components. For instance, extra "leaves" on the tree will be trimmed off when that tree is added to a phyloseq object.

**Example** - In the following example, we will define random example data tables in R, and then combine these into a phyloseq object. If you are able to get your data tables into R, then you can apply the following method to manually create phyloseq instances of your own data.


We'll create the example vanilla R tables using base R code. No packages required yet. 


```r
# Create a pretend OTU table that you read from a file, called otumat
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
##  [1,]   13   88   41   36   21   70   43   57   48    43
##  [2,]   19   56   79   54   67   22   50   69   56    61
##  [3,]   27   97   83   90   16   68   59   57   11    70
##  [4,]   11   82   90   85   97   30   16   46   80    23
##  [5,]   24   64   20   34   58   74   91   68   88    89
##  [6,]   93   35   41   79   87   13   99   87   28    74
##  [7,]    4   28   16   26   63   15   83   48   93    74
##  [8,]   73   67   87   49   43   58   48   37   15    63
##  [9,]   98   68   29   54   28   61   80   43   55    97
## [10,]   44    8   75   17   71   89   78   24   95    29
```


It needs sample names and OTU names, the index names of the your own matrix might already have this.


```r
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat
```

```
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       13      88      41      36      21      70      43      57
## OTU2       19      56      79      54      67      22      50      69
## OTU3       27      97      83      90      16      68      59      57
## OTU4       11      82      90      85      97      30      16      46
## OTU5       24      64      20      34      58      74      91      68
## OTU6       93      35      41      79      87      13      99      87
## OTU7        4      28      16      26      63      15      83      48
## OTU8       73      67      87      49      43      58      48      37
## OTU9       98      68      29      54      28      61      80      43
## OTU10      44       8      75      17      71      89      78      24
##       Sample9 Sample10
## OTU1       48       43
## OTU2       56       61
## OTU3       11       70
## OTU4       80       23
## OTU5       88       89
## OTU6       28       74
## OTU7       93       74
## OTU8       15       63
## OTU9       55       97
## OTU10      95       29
```


Now we need a pretend taxonomy table


```r
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
    "Species")
taxmat
```

```
##       Domain Phylum Class Order Family Genus Species
## OTU1  "j"    "m"    "w"   "e"   "a"    "i"   "e"    
## OTU2  "o"    "n"    "r"   "a"   "t"    "s"   "p"    
## OTU3  "x"    "o"    "d"   "x"   "s"    "z"   "x"    
## OTU4  "g"    "e"    "p"   "r"   "s"    "d"   "b"    
## OTU5  "i"    "h"    "l"   "j"   "p"    "b"   "g"    
## OTU6  "i"    "o"    "p"   "z"   "x"    "s"   "w"    
## OTU7  "b"    "d"    "c"   "m"   "j"    "e"   "v"    
## OTU8  "m"    "p"    "y"   "a"   "a"    "u"   "j"    
## OTU9  "y"    "d"    "v"   "o"   "z"    "g"   "m"    
## OTU10 "t"    "t"    "f"   "y"   "r"    "m"   "h"
```

```r
class(otumat)
```

```
## [1] "matrix"
```

```r
class(taxmat)
```

```
## [1] "matrix"
```


Note how these are just vanilla R matrices. Now let's tell phyloseq how to combine them into a phyloseq object.


```r
# In the previous lines, we didn't even need to have phyloseq loaded yet.
# Now we do.
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
```

```
## OTU Table:          [10 taxa and 10 samples]
##                      taxa are rows
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       13      88      41      36      21      70      43      57
## OTU2       19      56      79      54      67      22      50      69
## OTU3       27      97      83      90      16      68      59      57
## OTU4       11      82      90      85      97      30      16      46
## OTU5       24      64      20      34      58      74      91      68
## OTU6       93      35      41      79      87      13      99      87
## OTU7        4      28      16      26      63      15      83      48
## OTU8       73      67      87      49      43      58      48      37
## OTU9       98      68      29      54      28      61      80      43
## OTU10      44       8      75      17      71      89      78      24
##       Sample9 Sample10
## OTU1       48       43
## OTU2       56       61
## OTU3       11       70
## OTU4       80       23
## OTU5       88       89
## OTU6       28       74
## OTU7       93       74
## OTU8       15       63
## OTU9       55       97
## OTU10      95       29
```

```r
TAX
```

```
## Taxonomy Table:     [10 taxa by 7 taxonomic ranks]:
##       Domain Phylum Class Order Family Genus Species
## OTU1  "j"    "m"    "w"   "e"   "a"    "i"   "e"    
## OTU2  "o"    "n"    "r"   "a"   "t"    "s"   "p"    
## OTU3  "x"    "o"    "d"   "x"   "s"    "z"   "x"    
## OTU4  "g"    "e"    "p"   "r"   "s"    "d"   "b"    
## OTU5  "i"    "h"    "l"   "j"   "p"    "b"   "g"    
## OTU6  "i"    "o"    "p"   "z"   "x"    "s"   "w"    
## OTU7  "b"    "d"    "c"   "m"   "j"    "e"   "v"    
## OTU8  "m"    "p"    "y"   "a"   "a"    "u"   "j"    
## OTU9  "y"    "d"    "v"   "o"   "z"    "g"   "m"    
## OTU10 "t"    "t"    "f"   "y"   "r"    "m"   "h"
```

```r
physeq = phyloseq(OTU, TAX)
physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
```

```r
plot_bar(physeq, fill = "Family")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Let's add to this, pretending we also had other types of data available.

Create random sample data, and add that to the combined dataset. Make sure that the sample names match the `sample_names` of the `otu_table`.


```r
sampledata = sample_data(data.frame(Location = sample(LETTERS[1:4], size = nsamples(physeq), 
    replace = TRUE), Depth = sample(50:1000, size = nsamples(physeq), replace = TRUE), 
    row.names = sample_names(physeq), stringsAsFactors = FALSE))
sampledata
```

```
## Sample Data:        [10 samples by 2 sample variables]:
##          Location Depth
## Sample1         B   428
## Sample2         A   261
## Sample3         D   291
## Sample4         C   111
## Sample5         A   775
## Sample6         A    91
## Sample7         C   636
## Sample8         B   532
## Sample9         C   136
## Sample10        D   333
```


Now create a random phylogenetic tree with the ape package, and add it to your dataset. Make sure its tip labels match your `OTU_table`.


```r
library("ape")
random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Now let's combine these altogether. We can do this either by adding the new data components to the phyloseq object we already have by using `merge_phyloseq`, or we can use a fresh new call to `phyloseq` to build it again from scratch. The results should be identical, and we can check. You can always do either one with the help from accessor functions, and the choice is stylistic.

Merge new data with current phyloseq object:


```r
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]
```


Rebuild phyloseq data from scratch using all the simulated data components we just generated:


```r
physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
physeq2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]
```


Are they identical?


```r
identical(physeq1, physeq2)
```

```
## [1] TRUE
```


Let's build a couple tree plots with the new combined data.


```r
plot_tree(physeq1, color = "Location", label.tips = "taxa_names", ladderize = "left", 
    plot.margin = 0.3)
```

![plot of chunk treeplot](figure/treeplot1.png) 

```r
plot_tree(physeq1, color = "Depth", shape = "Location", label.tips = "taxa_names", 
    ladderize = "right", plot.margin = 0.3)
```

![plot of chunk treeplot](figure/treeplot2.png) 


Now how about some heatmaps.


```r
plot_heatmap(physeq1)
```

![plot of chunk heatmap-random](figure/heatmap-random1.png) 

```r
plot_heatmap(physeq1, taxa.label = "Phylum")
```

![plot of chunk heatmap-random](figure/heatmap-random2.png) 


As you can see, you gain access to the all the typical phyloseq tools, but without relying on any of the import wrappers.


---

### <a name="mgrast"></a> MG-RAST

A recent phyloseq issues tracker post discusses and demonstrates importing a `.biom` file exported by MG-RAST:

https://github.com/joey711/phyloseq/issues/272

The otherwise-recommended [import_biom](#import_biom) function does not work properly (for now) on this special variant of the BIOM-format. Or said another way, the [import_biom](#import_biom) function anticipates a different special variant of the BIOM-format the is generated by recent versions of QIIME. The [issue post about MG-RAST and phyloseq](https://github.com/joey711/phyloseq/issues/272) provides an example for [importing the data manually](#manual) using coercion functions and phyloseq constructors.


---

### <a name="microbio_me"></a> [microbio_me_qiime](http://joey711.github.io/phyloseq/download-microbio.me.html)

[microbio_me_qiime](http://joey711.github.io/phyloseq/download-microbio.me.html) is a function in phyloseq that interfaces with the QIIME data server:

http://www.microbio.me/qiime/index.psp

You will need to setup an account to browse the available data sets and their IDs. If you know a datasets ID already, or its assigned number, you can provide that as the sole argument to this function and it will download, unpack, and import the data into R, all in one command. Alternatively, if you have already downloaded the data from the QIIME server, and now have it locally on your hard drive, you can provide the local path to this tar-gz or zip file, and it will perform the unpacking and importing step for you. I'm finding this increasingly useful for creating demonstrations of methods and graphics, and can be a very effective way for you to provide fully reproducible analysis if your own data is already hosted on the [microbio.me](http://www.microbio.me/qiime/index.psp) server.

See the [microbio_me_qiime tutorial](http://joey711.github.io/phyloseq/download-microbio.me.html) for further details and examples.


---

### <a name="import_biom"></a>import_biom

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

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```r
plot_richness(myData, x = "BODY_SITE", color = "Description")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 

```r
plot_bar(myData, fill = "Genus")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) 

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

### <a name="import_qiime"></a> import_qiime

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

### <a name="import_mothur"></a> import_mothur


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
## otu_table()   OTU Table:         [ 591 taxa and 3 samples ]
## phy_tree()    Phylogenetic Tree: [ 591 tips and 590 internal nodes ]
```

```r
plot_tree(x, color = "samples")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-131.png) 

```r
SDF = data.frame(samples = sample_names(x), row.names = sample_names(x))
sample_data(x) = sample_data(SDF)
plot_richness(x)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-132.png) 


The class and data in the object returned by `import_mothur` depends on the  arguments. If the first three arguments are provided, then a phyloseq object should be returned containing both a tree and its associated OTU table. If only a list and group file are provided, then an "otu_table" object is returned. Similarly, if only a list and tree file are provided, then only a tree is returned ("phylo" class).

Returns just a tree

```r
x1 = import_mothur(mothlist, mothur_tree_file = mothtree, cutoff = "0.10")
x2 = import_mothur(mothlist, mothur_tree_file = mothtree, cutoff = "0.08")
plot(x1)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

Returns just an OTU table

```r
OTU = import_mothur(mothlist, mothgroup, cutoff = "0.08")
dim(OTU)
```

```
## [1] 591   3
```

```r
head(OTU)
```

```
## OTU Table:          [6 taxa and 3 samples]
##                      taxa are rows
##        B C D
## 9_6_14 2 0 0
## 9_1_14 1 0 0
## 9_1_15 1 0 0
## 9_1_16 1 0 0
## 9_1_18 1 0 0
## 9_1_19 1 0 0
```


Returns a list where each (outer) element represents an OTU, and is a vector of the sequencing reads that are clustered with that OTU.

```r
otulist = import_mothur(mothlist, cutoff = "0.08")
```

```
## Error: invalid subscript type 'list'
```

```r
length(otulist)
```

```
## Error: object 'otulist' not found
```

```r
ntaxa(OTU)
```

```
## [1] 591
```


Returns an error without a cutoff

```r
import_mothur(mothlist)
```

```
## Error: invalid subscript type 'list'
```


The list file is required. Import will fail with an error if it is not provided.

```r
import_mothur()
```

```
## Error: invalid subscript type 'list'
```



---

### <a name="import_pyrotagger"></a> import_pyrotagger

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

See [the tutorial on included example data in phyloseq](http://joey711.github.io/phyloseq/Example-Data) for more details.

The `data` command in the R language loads pre-imported datasets that are included in packages. For example, the "Global Patterns" dataset can be loaded into the R workspace with the following command.


```r
data(GlobalPatterns)
```




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

