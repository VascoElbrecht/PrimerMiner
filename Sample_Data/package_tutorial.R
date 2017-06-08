# Feel free to contact Vasco Elbrecht if you run into issues (twitter: @luckylionde). Enjoy!

# set the path to the PrimerMinder folder you just downloaded
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/PrimerMiner")

# install the PrimerMiner package icl dependencies
#install.packages(c("bold", "XML", "rentrez", "seqinr"), dependencies=T)
install.packages("PrimerMiner", repos = NULL, type="source", dependencies=T)

# load the package into R
library("PrimerMiner")

# quick guide inside R, further documentation in the Github Wiki
# https://github.com/VascoElbrecht/PrimerMiner/wiki
browseVignettes("PrimerMiner")


# Set path to sample data
setwd("Sample_Data")

# creating configuration file and batch downloading reads
batch_config("config.txt")

# batch download and process sequence data
batch_download("taxa_small.csv", "config.txt")

#You have to generate and manually check the alignment! Only extract the region amplified by the primers (+2 * the primer length on each side). We recommend Geneious for doing this!

# You find 5 sample alignments generated with Geneious in the folder "1 COI alignments (unprocessed)". We will use these in the next steps of this tutorial.

# get alignments to process
fastafiles <- list.files("1 COI alignments (unprocessed)", full.names=T)

# define name and folder to save files in
fastafiles_export <- paste("2 COI alignments (processed)", list.files("1 COI alignments (unprocessed)"), sep="/")

# process files! This function will remove gaps from the alignments and and apply selective trimming to the primer rgions
for (i in 1:length(fastafiles)){
selectivetrim(fastafiles[i], fastafiles_export[i], trimL=25, trimR=26, gaps=0.10, minsequL=100)
}


# The processed gap free alignemts can now be plotted with PrimerMiner or processed with third party primer development software
# Make sure to chekc that all alignemts have the same length (are aligned)
alignments <- list.files("2 COI alignments (processed)", full.name=T) # find files

pdf("plot_alignments.pdf", height=4, width=100)
plot_alignments(alignments, Order_names=gsub(".*./._(.*)_.*", "\\1", alignments))
dev.off()




# in silico primer evaluation
# With PrimerMiner v0.13 the scorring tables for missmatch possition and type are integrated in PrimerMiner. You can however generate your defult tables and use them (provide them as a csv, see sample data for the default tables!)


evaluate_primer("primer_scoring/01_Plecoptera_subset.fasta", "GGTCAACAAATCATAAAGATATTGG", 1, 25, save="save_evaluation_table_LCO.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 


evaluate_primer("primer_scoring/01_Plecoptera_subset.fasta", "ARYATDGTRATDGCHCCDGC", 585, 604, save="save_evaluation_table_BR1.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1", forward=F) 

# Evaluate primer paird against each others
primer_threshold("save_evaluation_table_LCO.csv", "save_evaluation_table_BR1.csv", threshold=120)





