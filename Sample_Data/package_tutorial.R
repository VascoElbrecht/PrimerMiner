# Feel free to contact Vasco Elbrecht if you run into issues (twitter: @luckylionde). Enjoy!


setwd("~/Desktop/PrimerMiner-master 4") # set the path to the PrimerMinder folder you just downloaded

# install the PrimerMiner package icl dependencies
install.packages("PrimerMiner", repos = NULL, type="source", dependencies=T)

# load the package into R
library("PrimerMiner")

setwd("Sample_Data")

# creating configuration file and batch downloading reads
batch_config("config.txt")

batch_download("taxa_small.csv", "config.txt")

##############
#Now you have to create alignments from the consensus OTU files! Please see manuals on github for more infos (not available jet!)
##############



# plot all alignments for manual primer finding!
alignments <- list.files("1 COi alignments", full.name=T) # find files

pdf("primer_plot_complete.pdf", height=6, width=100)
plot_alignments(alignments, Order_names=gsub(".*ts/.._(.*)_.*", "\\1", alignments))
dev.off()


# in silico primer evaluation

evaluate_primer("1 COi alignments/01_Plecoptera_folmer.fasta", "GCYCCHGAYATRGCHTTYCC", 218, 237, save="save_evaluation_table.csv", mm_position ="primer_scoring/Position_v1.csv", adjacent=2, mm_type="primer_scoring/Type_v1.csv") 








