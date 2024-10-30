Download_BOLD <- function(taxon, folder=NULL, setwd=NULL){

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}


cat(paste(Sys.time(), "- Downloading BOLD sequence data\n\n"), file= logfile, sep="", append=T)

folder_path <- ""
if(!is.null(folder)){dir.create(folder, showWarnings=F)
folder_path <- paste(folder, "/", sep="")
cat(paste("#Bold_data: Folder ",folder, "\n", sep=""), file= logfile, sep="", append=T)
} else {cat(paste("#Bold_data: ", "\n", sep=""), file= logfile, sep="", append=T)}
cat("Taxon\tSequences\tdownl_time\n", file= logfile, sep="", append=T)

# Use new BOLD library for boldsystems v5 
# library("bold")
if(!"BOLDconnectR" %in% loadedNamespaces()){
library("BOLDconnectR")
}
# Bold API key (log into your bold systems account to obtain)
if(exists("apikey")){
bold.apikey(apikey_bold) 
}

for (k in 1:length(taxon)){
time <- Sys.time() # get time
bold.data <- bold.public.search(taxonomy = taxon[k]) # fetch process IDs for specific taxons

data <- bold.fetch(get_by = "processid", identifiers = bold.data$processid)

# save tsv (without any processing)
if(save_bold_tsv){
write.table(data, file=paste(folder_path, taxon[k], "_BOLD.tsv"), quote=F, sep='\t', row.names=F)
}

# filter for only needed marker codes
data <- data[data$marker_code == marker_bold,]

# remove gabs from sequences
data$nuc <- gsub("-", "", data$nuc)
data <- data[nchar(data$nuc)>1,]


if(length(data)!=0){
cat(file=paste(folder_path, taxon[k], "_BOLD.fasta", sep="")) # delet old file
exp <- paste(">", data$processid, "___", data$order, "_", data$species, "\n", data$nuc, "\n", sep="")

cat(exp, file=paste(folder_path, taxon[k], "_BOLD.fasta", sep=""), append=T, sep="")
}
time <- Sys.time() - time
message(paste("Downloaded ", nrow(data)," sequences for ", taxon[k], " in ", format(time, digits=2), " from BOLD.", sep=""))
cat(paste(taxon[k],"\t", nrow(data), "\t", format(time, digits=2), "\n", sep=""), file= logfile, sep="", append=T)
}
cat("#Bold_data_end\n\n", file= logfile, sep="", append=T)
}


