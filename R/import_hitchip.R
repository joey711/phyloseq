#' Import HITChip output into phyloseq format
#'
#' @param data.dir Profiling script output directory for reading the data. 
#' @param method Probe summarization method ("rpa" or "sum")
#' @param detection.threshold Taxon absence/presence thresholds (typically 10^1.8 for HITChip)
#' @param verbose verbose
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples
#'   \dontrun{
#'     # Define example data folder with HITChip phylogenetic microarray data
#'     data.dir <- system.file("extdata/hitchip", package = "phyloseq")
#'     # Import the microarray data in phyloseq format
#'     physeq <- import_hitchip(data.dir)
#'   }
#'
#' @references 
#'   To cite the import_hitchip function, use
#'   "Leo Lahti 2015. import_hitchip function in the phyloseq R package."
#'   For the rpa method for probe summarization, cite:
#'   Lahti et al. NAR 41(10):e110, 2013.
#'   URL \url{http://nar.oxfordjournals.org/content/41/10/e110}
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
import_hitchip <- function(data.dir, method = "rpa", detection.threshold = 0, verbose = F) {

  # Read	     
  if ( verbose ) { message(paste("Reading Chip data from", data.dir)) }

  res <- list()

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  res[["probedata"]] <- tab

  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/taxonomy.filtered.tab", sep = "")    
  }
  taxonomy <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  res[["taxonomy"]] <- taxonomy
  
  # Read unfiltered taxonomy table
  f <- paste(data.dir, "/taxonomy.full.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/taxonomy.full.tab", sep = "")    
  }
  taxonomy.full <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  res[["taxonomy.full"]] <- taxonomy.full

  # Read sample metadata      
  f <- paste(data.dir, "/meta.tab", sep = "")
  if (file.exists(f)) {
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(tab) <- tab$sample
    meta <- tab
    res[["meta"]] <- meta
  }

  # -----------------------------------

  # Only pick probe-level data, taxonomy and metadata
  taxonomy <- res$taxonomy
  probedata <- res$probedata
  meta <- res$meta

  # Summarize probes into abundance table
  level <- "species"
  abu <- summarize_probedata(probedata = probedata,
      	 	             taxonomy = taxonomy,
      	 		     level = level,
			     method = method)

  # Convert the object into phyloseq format
  levels <- intersect(c("L0", "L1", "L2", "species"), colnames(taxonomy))
  taxonomy <- unique(taxonomy[, levels])
  rownames(taxonomy) <- taxonomy[, "species"]
  coms <- intersect(rownames(taxonomy), rownames(abu))
  abu <- abu[coms,]
  taxonomy <- taxonomy[coms,]

  pseq <- hitchip2physeq(t(abu), meta, taxonomy, detection.limit = detection.threshold)

  return(pseq) 
    
} 


#' Summarize phylogenetic microarray probe-level data from given input folder.
#' 
#' @param data.dir Data folder.
#' @param probedata probe-level data matrix
#' @param taxonomy probe taxonomy
#' @param level Summarization level
#' @param method Summarization method
#' @param probe.parameters Precalculater probe parameters. Optional.
#'
#' @return data matrix (taxa x samples)
#'
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
summarize_probedata <- function(data.dir = NULL, probedata = NULL, taxonomy = NULL, level, method, probe.parameters = NULL) {

  # If the data is not given as input, read it from the data directory
  if (method == "frpa" && is.null(probe.parameters)) {
     message("Loading pre-calculated RPA preprocessing parameters")
     probes <- unique(taxonomy[, "oligoID"])
     rpa.hitchip.species.probe.parameters <- list()
     load(system.file("extdata/probe.parameters.rda", package = "HITChipDB"))
     probe.parameters <- rpa.hitchip.species.probe.parameters
     # Ensure we use only those parameters that are in the filtered phylogeny
     for (bac in names(probe.parameters)) {
       probe.parameters[[bac]] <- probe.parameters[[bac]][intersect(names(probe.parameters[[bac]]), probes)]
     }
  }

  # Read probe-level data
  if (is.null(probedata)) {
    f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
    probedata <- tab
  }

  # Read taxonomy table
  if (is.null(taxonomy)) {
    f <- paste(data.dir, "/taxonomy.tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    # Convert into phyloseq taxonomyTable format
    taxonomy <- tax_table(as.matrix(tab))     
  }

  # Summarize probes through species level
  if (method %in% c("rpa", "frpa")) {
    otu <- summarize.rpa(taxonomy, level, probedata, verbose = TRUE, probe.parameters = probe.parameters)$abundance.table
  } else if (method == "sum") {
    otu <- summarize.sum(taxonomy, level, probedata, verbose = TRUE, downweight.ambiguous.probes = TRUE)$abundance.table
  }

  otu
    
}

#' Probeset summarization with SUM
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level taxonomic level for the summarization. 
#'   @param probedata preprocessed probes x samples data matrix in absolute domain
#'   @param verbose print intermediate messages
#'   @param downweight.ambiguous.probes Downweight probes with multiple targets
#'
#' Returns:
#'   @return List with two elements: abundance.table (summarized data matrix in absolute scale) and probe.parameters used in the calculations
#'
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
summarize.sum <- function (taxonomy, level, probedata, verbose = TRUE, downweight.ambiguous.probes = TRUE) {

  # Convert to log10 domain	      
  oligo.data <- probedata
  probe.parameters <- list()
 
  probesets <- retrieve.probesets(taxonomy, level = level)

  if (downweight.ambiguous.probes) {
    nPhylotypesPerOligo <- n.phylotypes.per.oligo(taxonomy, level) 
    probe.weights <- 1/nPhylotypesPerOligo
  } else {
    probe.weights <- rep(1, nrow(taxonomy))
    names(probe.weights) <- rownames(taxonomy)
  }

  # initialize
  summarized.matrix <- matrix(NA, nrow = length(probesets), 
  		       		  ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- names(probesets)
  colnames(summarized.matrix) <- colnames(oligo.data)

  for (set in names(probesets)) {

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- as.matrix(oligo.data[probes,])
    if (length(probes) == 1)  {
      vec <- as.vector(unlist(oligo.data[probes,]))
    } else {
      # Weight each probe by the inverse of the number of matching phylotypes
      # Then calculate sum -> less specific probes are downweighted
      # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      rownames(dat) <- probes
      colnames(dat) <- colnames(oligo.data)
      dat <- dat * probe.weights[rownames(dat)]
      vec <- colSums(dat, na.rm = T)               
    }

    summarized.matrix[set, ] <- vec

  }

  list(abundance.table = summarized.matrix, probe.parameters = probe.weights)
  
}


#' Description: Probeset summarization with RPA
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level taxonomic level for the summarization. 
#'   @param probedata preprocessed probes x samples data matrix in absolute domain
#'   @param verbose print intermediate messages
#'   @param probe.parameters Optional. If probe.parameters are given,
#'          the summarization is based on these and model parameters are not
#' 	    estimated. A list. One element for each probeset with the following probe vectors: 
#'	    affinities, variances
#' Returns:
#'   @return List with two elements: abundance.table (summarized data matrix in absolute scale) and probe.parameters (RPA probe level parameter estimates)
#'
#' @importFrom RPA d.update.fast 
#' @importFrom RPA rpa.fit
#'
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
summarize.rpa <- function (taxonomy, level, probedata, verbose = TRUE, probe.parameters = NULL) {

  # Convert to log10 domain	      
  oligo.data <- log10(probedata) 
  probeinfo <- list()
 
  probesets <- retrieve.probesets(taxonomy, level = level)
  # probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(taxonomy, level) 

  # initialize
  summarized.matrix <- matrix(NA, nrow = length(probesets), 
  		       		  ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- names(probesets)
  colnames(summarized.matrix) <- colnames(oligo.data)

  for (set in names(probesets)) {

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    if (length(probes) == 1)  {

      vec <- oligo.data[probes,]

    } else if (length(probe.parameters) > 0) {

      dat <- as.matrix(oligo.data[probes,])
      rownames(dat) <- probes
      colnames(dat) <- colnames(oligo.data)

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else {

      dat <- as.matrix(oligo.data[probes,])
      rownames(dat) <- probes
      colnames(dat) <- colnames(oligo.data)

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      alpha <- 1 + 0.1*ncol(dat)/2
      beta <- 1 + 0.1*ncol(dat)*nPhylotypesPerOligo[probes]^2
      res <- rpa.fit(dat, alpha = alpha, beta = beta)
      vec <- res$mu
      probeinfo[[set]] <- res$tau2

    }
      
    summarized.matrix[set, ] <- as.vector(unlist(vec))

  }

  if (!is.null(probe.parameters)) {
    probeinfo <- probe.parameters
  }

  # Return the data in absolute scale					
  summarized.matrix <- 10^summarized.matrix

  list(abundance.table = summarized.matrix, probeinfo = probeinfo)
  
}


#' hitchip2physeq
#'
#' Convert HITChip data into phyloseq format
#'
#' @param otu Sample x OTU absolute HITChip signal
#' @param meta Sample x features metadata data.frame
#' @param taxonomy OTU x Taxonomy data.frame (HITChip taxonomy used by default)
#' @param detection.limit HITChip signal detection limit (absence / presence)
#' @return phyloseq object
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
hitchip2physeq <- function (otu, meta, taxonomy = NULL, detection.limit = 10^1.8) {

  # OTU x Sample matrix: absolute 'read counts'
  x <- t(otu) - detection.limit # HITChip detection limit
  x[x < 0] <- 0
  x <- 1 + x

  # -------------------------

  # Discretize to get 'counts'
  otumat <- round(x)
  OTU <- otu_table(otumat, taxa_are_rows = TRUE)

  # Create phyloseq object
  pseq <- phyloseq(OTU)

  # --------------------------

  # Construct taxonomy table
  # If nrow(otumat) then it is probe-level data and no taxonomy should be given
  # for that by default
  if (is.null(taxonomy) && nrow(otumat) < 3000) {

    ph <- GetPhylogeny("HITChip")
    ph <- unique(ph[, c("L1", "L2", "species")])
    colnames(ph) <- c("Phylum", "Genus", "Phylotype")
    taxonomy <- ph
    input.level <- colnames(ph)[[which.max(apply(ph, 2, function (x) {sum(rownames(otumat) %in% x)}))]]
    if (input.level == "Genus") {
      taxonomy <- unique(taxonomy[, c("Phylum", "Genus")])
    } else if (input.level == "Phylum") {
      taxonomy <- as.data.frame(Phylum = unique(taxonomy[, c("Phylum")]), ncol = 1)
      	       
    }
    rownames(taxonomy) <- as.character(taxonomy[[input.level]])
  }

  if (!all(rownames(otumat) %in% rownames(taxonomy)) && nrow(otumat) < 3000) {
      warning(paste("Some OTUs are missing from the taxonomy tree!", paste(setdiff(rownames(otumat), rownames(taxonomy)), collapse = " / ")))
      # Common probes or OTUs
      coms <- intersect(rownames(otumat), rownames(taxonomy))
      # Only keep probes that have taxonomy information
      otumat <- otumat[coms, ]
      taxonomy <- taxonomy[coms, ]
   }

   if (!is.null(taxonomy) || nrow(otumat) < 3000) {
     TAX <- tax_table(as.matrix(taxonomy[rownames(otumat), ]))

     # Combine OTU and Taxon matrix into Phyloseq object
     pseq <- merge_phyloseq(pseq, TAX)
  }
  
  # -------------------------

  if (!is.null(meta)) {
  
    # Metadata
    rownames(meta) <- as.character(meta$sample)
    sampledata <- sample_data(meta[colnames(otumat),])

    pseq <- merge_phyloseq(pseq, sampledata)

    # Harmonize the fields?
    # pseq@sam_data <- harmonize_fields(pseq@sam_data)
    
  }
  
  # --------------------------

  # We could also add phylotree between OTUs
  # source("tree.R")
  # pseq <- merge_phyloseq(pseq, tree2)
  # pseq <- merge_phyloseq(pseq, random_tree)

  # --------------------------

  pseq
 
}


#' retrieve.probesets
#' 
#' List probes for each probeset
#'
#' @param tax.table data.frame with oligo - phylotype 
#' 	  		 mapping information
#' @param level phylotype level for probesets
#' @param name specify phylotypes to check (optional)
#'
#' @return A list. Probes for each phylotype.
#'
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
retrieve.probesets <- function(tax.table, level = "species", name = NULL) {

    tax.table <- as.data.frame(tax.table)

    # If name not given, pick all
    if (is.null(name)) {
        name <- unique(as.character(tax.table[[level]]))
    }
    
    phylo <- tax.table[tax.table[[level]] %in% name, ]
    
    if (is.factor(phylo[[level]])) {
        phylo[[level]] <- droplevels(phylo[[level]])
    }
    
    phylo.list <- split(phylo, phylo[[level]])
    probesets <- lapply(phylo.list, function(x) {
        as.character(unique(x$oligoID))
    })
    names(probesets) <- names(phylo.list)
    
    probesets
    
} 



#' Description: Check number of matching phylotypes for each probe
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level phylotype level
#'
#' Returns:
#'   @return number of matching phylotypes for each probe
#'
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal
n.phylotypes.per.oligo <- function (taxonomy, level) {
  sapply(split(as.vector(taxonomy[, level]), as.vector(taxonomy[, "oligoID"])), function(x) length(unique(x)))
}

