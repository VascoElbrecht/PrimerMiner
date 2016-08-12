# 150223


download_stats <- function(table, save="Download_Stats_table.csv", cluster="Download_Stats_cluster.csv") {

if(is.data.frame(table)){} else {table <- read.csv(table, sep= Taxon_sep, stringsAsFactors=F)}


folder <- which(table[,1]!="")
folder <- c(folder, nrow(table)+1)


if (cluster!=""){
Vstat <- data.frame(Order=table[table[,1]!="",1], Sequ=NA, Derep=NA, Cluster=NA)}

i <- 1
for (i in 1:(length(folder)-1)){
subFolder <- table[folder[i],1]
subStart <- folder[i]+1
subEnd <- folder[i+1]-1

if (subStart > subEnd) {subStart <- subEnd}

log <- readLines(paste(subFolder, "/log.txt", sep=""))

if(download_bold){
BOLD <- grep("#Bold_data", log)[1:2]
BOLD <- log[(BOLD[1]+2):(BOLD[2]-1)]
BOLD <- as.numeric(sub(".*\t(.*)\t.*", "\\1", BOLD))
table[subStart:subEnd,3] <- BOLD
if (is.na(table[folder[i],3])) {table[folder[i],3] <- sum(BOLD)}
} else {table[folder[i],3] <- NA}


if(download_GB){
GB <- grep("#GB_data", log)[1:2]
GB <- log[(GB[1]+2):(GB[2]-1)]
GB <- as.numeric(sub(".*\t(.*)\t.*", "\\1", GB))
table[subStart:subEnd,4] <- GB
if (is.na(table[folder[i],4])) {table[folder[i],4] <- sum(GB)}
} else {table[folder[i],4] <- NA}



if(download_mt){
mito <- grep("#mito_data", log)[1:2]
mito <- log[(mito[1]+2):(mito[2]-1)]
mito <- as.numeric(sub(".*\t(.*)\t.*", "\\1", mito))
table[subStart:subEnd,5] <- mito
if (is.na(table[folder[i],5])) {table[folder[i],5] <- sum(mito)}
} else {table[folder[i],5] <- NA}


# cluster stats
if (cluster!=""){
imp <- sum(table[folder[i],3], table[folder[i],4], table[folder[i],5], na.rm=T)
#as.numeric(sub(".* (.*)$", "\\1", log[grep("Number of imput sequences", log)[1]]))
derep <- as.numeric(sub(".* (.*)$", "\\1", log[grep("Dereplicated: ", log)[1]]))
clust <- as.numeric(sub(".* (.*)$", "\\1", log[grep("Cluster: ", log)[1]]))

Vstat[i,-1] <- c(imp, derep, clust)
}

}
names(table)[3:5] <- c("BOLD", "GB", "Mito")
write.table(table, save, sep="\t", row.names=F)
message(paste("Download stats written in table ", save, sep=""))


# cluster stats
if (cluster!=""){
write.table(Vstat, file=cluster, sep="\t", row.names=F)
message(paste("Cluster stats written in table ", cluster, sep=""))

Vstat[i,-1] <- c(imp, derep, cluster)
}


}


