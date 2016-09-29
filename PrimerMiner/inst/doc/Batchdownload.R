## ---- eval = F-----------------------------------------------------------
#  setwd("your/path/to/PrimerMiner/Sample_Data")

## ---- eval = F-----------------------------------------------------------
#  batch_config("config.txt")

## ---- eval = F-----------------------------------------------------------
#  batch_download("taxa_small.csv", "config.txt")

## ---- eval = F-----------------------------------------------------------
#  fastafiles <- list.files("1 COI alignments (unprocessed)", full.names=T)
#  fastafiles_export <- paste("2 COI alignments (processed)", list.files("1 COI alignments (unprocessed)"), sep="/")

## ---- eval = F-----------------------------------------------------------
#  for (i in 1:length(fastafiles)){
#  selectivetrim(fastafiles[i], fastafiles_export[i], trimL=25, trimR=26, gaps=0.10, minsequL=100)
#  }

## ---- eval = F-----------------------------------------------------------
#  alignments <- list.files("2 COI alignments (processed)", full.name=T)
#  
#  pdf("plot_alignments.pdf", height=4, width=100)
#  plot_alignments(alignments, Order_names=gsub(".*./._(.*)_.*", "\\1", alignments))
#  dev.off()

