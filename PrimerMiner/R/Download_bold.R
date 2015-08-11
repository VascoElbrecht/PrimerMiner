Download_BOLD <- function(taxon, folder=NULL, setwd=NULL){

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}


cat(paste(Sys.time(), "- Downloading BOLD sequence data\n\n"), file= logfile, sep="", append=T)

folder_path <- ""
if(!is.null(folder)){dir.create(folder, showWarnings=F)
folder_path <- paste(folder, "/", sep="")
cat(paste("#Bold_data: Folder ",folder, "\n", sep=""), file= logfile, sep="", append=T)
} else {cat(paste("#Bold_data: ", "\n", sep=""), file= logfile, sep="", append=T)}
cat("Taxon\tSequences\tdownl_time\n", file= logfile, sep="", append=T)

library("bold")

for (k in 1:length(taxon)){
time <- Sys.time() # get time
data <- bold_seq(taxon=taxon[k])

if(length(data)!=0){
cat(file=paste(folder_path, taxon[k], "_BOLD.fasta", sep="")) # delet old file
for (i in 1:length(data)){
exp <- paste(">", data[i][[1]][1], "_", data[i][[1]][2], "\n", gsub("-", "", data[i][[1]][4]), "\n", sep="")
cat(exp, file=paste(folder_path, taxon[k], "_BOLD.fasta", sep=""), append=T)
}
}
time <- Sys.time() - time
print(paste("Downloaded ", length(data)," sequences for ", taxon[k], " in ", format(time, digits=2), " from BOLD.", sep=""))
cat(paste(taxon[k],"\t", length(data), "\t", format(time, digits=2), "\n", sep=""), file= logfile, sep="", append=T)
}
cat("#Bold_data_end\n\n", file= logfile, sep="", append=T)
}


