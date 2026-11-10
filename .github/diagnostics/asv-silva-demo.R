# Entirely synthetic offline integration demo. Not a biological benchmark.
# Run from the repository root; use a fresh output directory.
args <- commandArgs(trailingOnly=TRUE)
output <- if (length(args)) args[[1]] else tempfile("asv-silva-demo-")
if (dir.exists(output)) stop("Choose a new output directory.")
dir.create(output,recursive=TRUE)
pkgload::load_all(".",quiet=TRUE)
set.seed(942)
base <- sample(c("A","C","G","T"),1200,replace=TRUE)
references <- vapply(seq_len(4),function(i) {
  sequence <- base
  positions <- seq(i*7,1200,by=83)
  sequence[positions] <- sample(c("A","C","G","T"),length(positions),TRUE)
  paste(sequence,collapse="")
},character(1))
accessions <- paste0("synthetic",seq_along(references))
ids <- paste0(accessions,".1.1200")
names(references) <- paste(ids,"Bacteria;Synthetic offline fixture")
fasta <- file.path(output,"synthetic-reference.fasta")
quality_file <- file.path(output,"synthetic-quality.tsv")
Biostrings::writeXStringSet(Biostrings::DNAStringSet(references),fasta)
quality <- data.frame("Primary Accession"=accessions,Start=1,Stop=1200,
                        "Sequence Quality"=100,"Alignment Quality"=100,
                        "% Ambiguities"=0,"% Homopolymers"=0,"% Vector Contamination"=0,
                        check.names=FALSE)
write.table(quality,quality_file,sep="\t",row.names=FALSE,quote=FALSE)
queries <- c("01"=substr(references[[1]],300,549),
               "asv2"=substr(references[[3]],350,599))
counts <- matrix(c(20,2,7,18),2,dimnames=list(names(queries),c("01","sample2")))
input <- phyloseq(otu_table(counts,TRUE),Biostrings::DNAStringSet(queries))
original <- input
selected <- select_asvs(input,ids=names(queries))
# Version here describes the fixture's metadata schema, not real SILVA content.
raw <- silva_reference("144",fasta=fasta,quality=quality_file)
filtered <- filter_silva_reference(raw,file.path(output,"qc"),chunk_size=2)
hits <- match_asvs_to_reference(selected,filtered,max_mismatches=8,chunk_size=2)
result <- infer_asv_tree(hits,model="GTR",bootstrap=0,seed=42)
stopifnot(identical(input,original),filtered$summary$retained==4,
          setequal(result$asv_tree$tip.label,names(queries)))
saveRDS(result,file.path(output,"phylogeny.rds"))
Biostrings::writeXStringSet(result$alignment,file.path(output,"alignment.fasta"))
write.table(hits$matches,file.path(output,"candidates.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(hits$summary,file.path(output,"query-summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Newick uses safe internal labels; tip-data.tsv supplies the original-ID map.
ape::write.tree(result$fit$tree,file=file.path(output,"context-internal-ids.nwk"))
write.table(result$tip_data,file.path(output,"tip-data.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
origin <- result$tip_data$Origin[match(result$tree$tip.label,result$tip_data$TipLabel)]
grDevices::png(file.path(output,"context-tree.png"),width=1400,height=850,res=140)
par(mar=c(4,1,4,1))
plot(result$tree,tip.color=ifelse(origin=="ASV","#0072B2","#555555"),cex=0.9,
     font=1,y.lim=c(0,7),no.margin=FALSE)
title(paste("Synthetic 16S reference context and observed ASV fragments",
            "Joint GTR ML tree; unrooted; bootstrap support not assessed",sep="\n"),cex.main=0.9)
ape::add.scale.bar(y=0.25,cex=0.8)
legend("bottomright",legend=c("Observed ASV (250 nt)","Synthetic reference (1200 nt)"),
       text.col=c("#0072B2","#555555"),bty="n",cex=0.8)
invisible(grDevices::dev.off())
writeLines(c("Synthetic offline example; not SILVA data or a biological benchmark.",
               capture.output(sessionInfo())),file.path(output,"session-info.txt"))
cat("Offline demonstration saved to:",normalizePath(output),"\n")
