Download_GB <- function(taxon, folder=NULL, marker=c("COi", "CO1", "COXi", "COX1"), maxlength=2000, custom_query=NULL, setwd=NULL){

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}


cat(paste(Sys.time(), "- Downloading GenBank sequence data\n\n"), file=logfile, sep="", append=T)

if (is.null(custom_query)){cat(paste("Search query: REPLACE_WITH_TAXA", "[Organism] AND (", paste(c(marker), collapse=" OR "), ") AND 1:",maxlength ,"[Sequence Length]\n\n", sep=""), file=logfile, sep="", append=T)} else {cat(paste("Search query: REPLACE_WITH_TAXA", custom_query, "\n\n", sep=""), file=logfile, sep="", append=T)}

folder_path <- ""
if(!is.null(folder)){dir.create(folder, showWarnings=F)
folder_path <- paste(folder, "/", sep="")
cat(paste("#GB_data: Folder ",folder, "\n", sep=""), file=logfile, sep="", append=T)
} else {cat(paste("#GB_data: ", "\n", sep=""), file=logfile, sep="", append=T)}
cat("Taxon\tSequences\tdownl_time\n", file=logfile, sep="", append=T)

for (k in 1:length(taxon)){
time <- Sys.time() # get time


# download IDs
if (is.null(custom_query)){
searchQ <- paste(taxon[k],"[Organism] AND (", paste(c(marker), collapse=" OR "), ") AND 1:",maxlength ,"[Sequence Length]", sep="")
} else {searchQ <- paste(taxon, custom_query, sep="")}

search_results <- entrez_search(db="nuccore", term=searchQ, retmax=9999999)


if(length(search_results$ids)!=0){
	
cat("", file=paste(folder_path, taxon[k], "_GB.fasta", sep=""), sep="", append=F) # overwrite old file!


i <- 1
while (!is.na(search_results$ids[i])){
temp <- search_results$ids[i:(i+499)]
temp <- temp[!is.na(temp)]
downloaded_sequ <- entrez_fetch(db="nuccore", id=temp, rettype="fasta")
if (downloaded_sequ[1]!="resource temporarily unavailable (4)."){
cat(downloaded_sequ, file=paste(folder_path, taxon[k], "_GB.fasta", sep=""), sep="", append=T)
i <- i + 500} # only write in file if data downloaded!
Sys.sleep(0.5)
}

meep <- read.fasta(paste(folder_path, taxon[k], "_GB.fasta", sep=""))
if (length(meep)!=length(search_results$ids)){
paste("WARNING: Something went wrong with the download. Numer of files in GB does not match number of downloaded files!")
cat(paste("\nWARNING: Something went wrong with the download. Numer of files in GB does not match number of downloaded files!\n\n"), file=logfile, sep="", append=T)}

}
time <- Sys.time() - time
print(paste("Downloaded ", length(search_results$ids)," sequences for ", taxon[k], " in ", format(time, digits=2), " from NCBI.", sep=""))
cat(paste(taxon[k],"\t", length(search_results$ids), "\t", format(time, digits=2), "\n", sep=""), file=logfile, sep="", append=T)
}

cat("#GB_data_end\n\n", file=logfile, sep="", append=T)
}

