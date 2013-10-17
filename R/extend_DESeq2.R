phyloseq2deseq2 = function(physeq, ...){
  require("DESeq2")
  require("phyloseq")
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  countData = countData + 1L
  # Need to add check here for missing sample_data
  colData = data.frame(sample_data(physeq))
  # Re-order the levels so the NULL set is first.
  colData$postfix <- relevel(colData$postfix, levels(colData$postfix)[2])
  
  # Need to add taxonomyTable support...
  
  # Create the DESeq data set, dds.
  dds <- DESeqDataSetFromMatrix(countData, colData, ...)
  return(dds)
}