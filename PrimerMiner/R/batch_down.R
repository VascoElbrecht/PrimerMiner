batch_download <- function(table, config){


source("config.txt")
if(is.data.frame(table)){} else {table <- read.csv(table, sep= Taxon_sep, stringsAsFactors=F)}


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
print(paste("Data for *", subFolder, "* was already downloaded and clustered on ", time[1], " and will thus be skipped. Turn Skip_if_complete to F, if you like data to be redownloaded and reclustered or delete the file done.txt in the folder ", subFolder,  sep=""))
print("")
print("-------------------")
print("")
}
}

if(download_and_cluster){

# download data!

if (Download){
if (download_bold){Download_BOLD(taxa, folder= subFolderPath, setwd=subFolder)}
if (download_GB){Download_GB(taxa, folder= subFolderPath, marker=Marker, maxlength= maxlength_GB, custom_query= custom_query_GB, setwd=subFolder)}
if (download_mt){Download_mito(taxa, folder= subFolderPath, minlength= minlength_mt, maxlength= maxlength_mt, custom_query = custom_query_mt, setwd=subFolder)
	if (length(list.files(subFolderPath, pattern="mito.gb$"))>0){
		print("")
		Mito_GB2fasta(subFolderPath, marker=Marker, add= add_mt, rm_dup= rm_dup, no_marker= no_marker, setwd=subFolder)}
}
}
print("")
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


Merge_fasta(all_file_TF, save=all_fasta , clip_left=0, clip_right=0, setwd=subFolder)


print("")


all_fasta <- paste(subFolder, "_all.fasta", sep="")

Clustering(all_fasta, vsearchpath= vsearchpath, id=id, threshold=threshold, setwd=subFolder)
print("")
print("-------------------")
print("")


cat(file=paste(subFolder, "/", "done.txt", sep=""), paste(Sys.time()))

}
}# download and cluster or not!

} # folder end

#print("starting with stats")
# write summary statistics
if(summstats){
download_stats(table)
}

print("")
print("We are all done here. Please use a porgram like Geneious to create sequence  alignments. reffer to the manual on GitHub for further instructions!")
print("https://github.com/VascoElbrecht/PrimerMiner")
} # function end




