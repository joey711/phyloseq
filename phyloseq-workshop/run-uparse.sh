#!/bin/sh

#  run-uparse.sh
#  
#
#  Created by Paul Joseph McMurdie II on 5/22/14.
#  Modified on 2014-07-08
#
# Feel free to also check out Robert Edgar's own suggested "pipeline"
# of BASH commands:
# http://drive5.com/usearch/manual/uparse_cmds.html
#
################################################################################
#
# run-uparse.sh
# You need to change these names for your specific system/run
# Also modify file output names as needed/desired. Make sure you
# use find/replace while modifying a file name, as an output at one step
# might be an input at a later step...
#
################################################################################
# Define machine-specific path variables
# You should change these according to your machine.

# Set a variable to give a short name for the USEARCH binary
# Don't forget that it needs to be executable
# e.g. run `chmod +x usearch7.0.1090_i86osx32` when you first set it up.
u=./usearch7.0.1090_i86osx32

# Location of auxiliary python scripts
# http://drive5.com/python/python_scripts.tar.gz
# 
p=python_scripts

# Define gold reference database path for chimera checking
# http://drive5.com/uchime/uchime_download.html
# http://drive5.com/uchime/gold.fa
gold=gold.fa

# Set a variable name for original input fasta (.fna) file
i=qiime_overview_tutorial/seqs-trim-210.fna

# Define path to reference database (e.g. greengenes)
gg=greengenes/13_8_97_otus.fasta

################################################################################
#
# Dereplication, chimera-filtering, and de-novo clustering
#
################################################################################

echo "Dereplication..."
$u -derep_fulllength $i -output dereplicated.fna -sizeout

echo "Abundance sort and discard singletons..."
$u -sortbysize dereplicated.fna -output dereplicated_sorted.fna -minsize 2

echo "De-Novo OTU clustering (and denovo chimera removal)"
$u -cluster_otus dereplicated_sorted.fna -otus denovo_no_chimera.fna

echo "Chimera filtering using reference database"
$u -uchime_ref denovo_no_chimera.fna -db $gold -strand plus -nonchimeras uparse_OTUs.fna

echo "Labeling OTUs and matching against the gg database for closed and open results..."
# Label OTU sequences
python $p/fasta_number.py uparse_OTUs.fna OTU_ > uparse_OTUs_Labelled.fna

# Label OTUs accordingly if a match against the greengenes OTUs database
# http://greengenes.secondgenome.com/downloads/database/13_5
$u -usearch_global uparse_OTUs_Labelled.fna -db $gg -strand both -id 0.97 -dbmatched ggmatch_OTUs.fna -notmatched novel_OTUs.fna

# Combined the gg OTUs that matched, and the novel OTUs that didn't
cat ggmatch_OTUs.fna novel_OTUs.fna > open_OTUs.fna

################################################################################
#
# De Novo
#
################################################################################
echo "Map original input reads (including singletons) back to denovo OTUs..."
$u -usearch_global $i -db uparse_OTUs_Labelled.fna -strand plus -id 0.97 -uc denovo_map.uc

################################################################################
#
# Open-Reference
#
################################################################################
echo "Map original input reads (including singletons) back to open OTUs..."
$u -usearch_global $i -db open_OTUs.fna -strand plus -id 0.97 -uc open_map.uc

################################################################################
#
# Closed-Reference
#
################################################################################
echo "Map our original reads (including singletons) against..."
$u -usearch_global $i -db ggmatch_OTUs.fna -strand both -id 0.97 -uc closed_map.uc

echo "DONE! :-)"
