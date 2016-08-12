Download_mito <- function(taxon, folder=NULL, minlength=2001, maxlength=80000, custom_query=NULL, setwd=NULL){

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}


cat(paste(Sys.time(), " - Downloading Miochondrial Genomes from GenBank\n\n"), file=logfile, sep="", append=T)

if (is.null(custom_query)){cat(paste("Search query: REPLACE_WITH_TAXA", "[Organism] AND mitochondrion[filter] AND genome AND ", minlength, ":", maxlength ,"[Sequence Length]\n\n", sep=""), file=logfile, sep="", append=T)} else {cat(paste("Search query: REPLACE_WITH_TAXA", custom_query, "\n\n", sep=""), file=logfile, sep="", append=T)}

folder_path <- ""
if(!is.null(folder)){dir.create(folder, showWarnings=F)
folder_path <- paste(folder, "/", sep="")
cat(paste("#mito_data: Folder ",folder, "\n", sep=""), file=logfile, sep="", append=T)
} else {cat(paste("#mito_data: ", "\n", sep=""), file=logfile, sep="", append=T)}
cat("Taxon\tSequences\tdownl_time\n", file=logfile, sep="", append=T)

for (k in 1:length(taxon)){
time <- Sys.time() # get time

# download IDs
if (is.null(custom_query)){
searchQ <- paste(taxon[k],"[Organism] AND mitochondrion[filter] AND genome AND ", minlength, ":", maxlength ,"[Sequence Length]", sep="")
} else {searchQ <- paste(taxon, custom_query, sep="")}

search_results <- entrez_search(db="nuccore", term=searchQ, retmax=9999999)

# save genbank file!
if(length(search_results$ids)!=0){
cat(file=paste(folder_path, taxon[k], "_mito.gb", sep="")) # overwrite old files


for (i in 1:length(search_results$ids)){
download.file(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", search_results$ids[i], "&rettype=gb&retmode=text", sep=""), destfile=paste(folder_path, taxon[k], "_mito.gb", sep=""), mode="a", quiet=T)}

}

time <- Sys.time() - time
print(paste("Downloaded ????????", length(search_results$ids)," mitogenomes for ", taxon[k], " in ", format(time, digits=2), ".", sep=""), file=logfile, sep="", append=T)
cat(paste(taxon[k],"\t", length(search_results$ids), "\t", format(time, digits=2), "\n", sep=""), file=logfile, sep="", append=T)

}
cat(paste("#mito_data_end", "\n\n", sep=""), file=logfile, sep="", append=T)
#return(length(search_results$ids))
}

