# 150324 - test file!
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/PrimerMiner/") # set the path to the PrimerMinder folder you just downloaded

# install packages needed for PrimerMiner
install.packages(c("bold", "XML", "rentrez", "seqinr"), dependencies=T)

# install the PrimerMiner package
install.packages("PrimerMiner", repos = NULL, type="source")

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
alignments <- list.files("Sample_Data/1 COi alignments", full.name=T) # find files

pdf("Sample_Data/primer_plot_complete.pdf", height=6, width=100)
plot_alignments(alignments, Order_names=gsub(".*ts/.._(.*)_.*", "\\1", alignments))
dev.off()


