batch_download <- function(table, config){


source(config)
bold.apikey(apikey_bold) 

if(is.data.frame(table)){} else {table <- read.csv(table, sep= Taxon_sep, stringsAsFactors=F)}

table[2][is.na(table[2])] <- "" # replace NAs if only orders are given

folder <- which(table[,1]!="")
folder <- c(folder, nrow(table)+1)

table[folder,]

for (i in 1:(length(folder)-1)){
subFolder <- table[folder[i],1]
subStart <- folder[i]+1
subEnd <- folder[i+1]-1

dir.create(subFolder, showWarnings=F)
#setwd(subFolder) # set up folder and paste files!

subFolderPath <- paste(subFolder, "/", subFolder, sep="") # set op filepath



# define range!
if (table[subEnd,2]==""){taxa <- subFolder} else {
taxa <- table[subStart:subEnd, 2]}


# check if downloading and clustering is already done!
download_and_cluster <- T # tipicly data has not been processed

if(Skip_if_complete){
done <- file.exists(paste(subFolder, "/", "done.txt", sep=""))

if(done){
download_and_cluster <- F # DON'T REDOWNLOAD DATA
time <- readLines(paste(subFolder, "/", "done.txt", sep=""), warn=F)
message(paste("Data for *", subFolder, "* was already downloaded and clustered on ", time[1], " and will thus be skipped. Turn Skip_if_complete to F, if you like data to be redownloaded and reclustered or delete the file done.txt in the folder ", subFolder,  sep=""))
message(" ")
message("-------------------")
message(" ")
}
}

if(download_and_cluster){

# download data!

if (Download){
message("Starting to downloand data for ", subFolder, ": \"", paste(taxa, collapse="\", \""), "\"")
if (download_bold){Download_BOLD(taxa, folder= subFolderPath, setwd=subFolder)}
if (download_GB){Download_GB(taxa, folder= subFolderPath, marker=Marker, maxlength= maxlength_GB, custom_query= custom_query_GB, setwd=subFolder, GB_subset=GB_Subset)}
if (download_mt){Download_mito(taxa, folder= subFolderPath, minlength= minlength_mt, maxlength= maxlength_mt, custom_query = custom_query_mt, setwd=subFolder, MT_Subset= mt_subset)
	if (length(list.files(subFolderPath, pattern="mito.gb$"))>0){
		message(" ")
		Mito_GB2fasta(subFolderPath, marker=Marker, add= add_mt, rm_dup= rm_dup, no_marker= no_marker, setwd=subFolder)}
}
}
message(" ")
# merge files!
if (Merge_and_Cluster_data){


all_file_TF <- c()

if (merge_bold){
BOLD_fasta <- paste(subFolderPath, "_Bold.fasta", sep="")

if(length(list.files(subFolderPath, pattern="BOLD\\.fasta$"))>0){
Merge_fasta(files=subFolderPath, save= BOLD_fasta, clip_left= clipping_left_bold, clip_right= clipping_rigth_bold, pattern="BOLD\\.fasta$", setwd=subFolder)
all_file_TF <- c(all_file_TF, BOLD_fasta)}
}

if (merge_GB){
if(length(list.files(subFolderPath, pattern="GB\\.fasta$"))>0){
GB_fasta <- paste(subFolderPath, "_GB.fasta", sep="")
Merge_fasta(subFolderPath, save= GB_fasta, clip_left= clipping_left_GB, clip_right= clipping_rigth_GB, pattern="GB\\.fasta$", setwd=subFolder)
all_file_TF <- c(all_file_TF, GB_fasta)}
}

if (merge_mt) {
if(length(list.files(subFolderPath, pattern="[mito]\\.fasta$"))>0){
mito_fasta <- paste(subFolderPath, "_mito.fasta", sep="")
Merge_fasta(subFolderPath, save= mito_fasta, clip_left= clipping_left_mt, clip_right= clipping_rigth_mt, pattern="[mito]\\.fasta$", setwd=subFolder)
all_file_TF <- c(all_file_TF, mito_fasta)}
}

# all files together
all_fasta <- paste(subFolder, "/", subFolder, "_all.fasta", sep="")

if(is.null(all_file_TF)){
glumanda <- paste("\nWARNING: For the group ", subFolder, " no sequences were obtained from the given reference databases. review the taxon spelling or search for a broader group / aktivate downloading on all databases.\n\n")
cat(file=paste(subFolder, "/log.txt", sep=""), glumanda, append=T)
message(glumanda)
} else {

Merge_fasta(all_file_TF, save=all_fasta , clip_left=0, clip_right=0, setwd=subFolder)


message(" ")


all_fasta <- paste(subFolder, "_all.fasta", sep="")

Clustering(all_fasta, vsearchpath= vsearchpath, id=id, threshold=threshold, setwd=subFolder)
message(" ")
} # empy file else done
message("-------------------")
message(" ")


cat(file=paste(subFolder, "/", "done.txt", sep=""), paste(Sys.time()))

}
}# download and cluster or not!

} # folder end

#print("starting with stats")
# write summary statistics
if(summstats){
download_stats(table)
}

message(" ")
message("We are all done here. Please use a porgram like Geneious to create sequence  alignments. reffer to the wiki on GitHub for further instructions!")
message("https://github.com/VascoElbrecht/PrimerMiner/wiki")
} # function end




