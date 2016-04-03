# Vsearch


Clustering <- function(file, vsearchpath="integrated", id=0.97, cmd="", threshold="Majority", setwd=NULL){


if (vsearchpath=="integrated"){
if (operating_system=="MacOSX"){
vsearchpath <- paste(system.file(package="PrimerMiner"), "/vsearch-1.10.2_osx_x86_64", sep="")
Sys.chmod(vsearchpath, mode = "0777", use_umask = TRUE)} else {
vsearchpath <- paste(system.file(package="PrimerMiner"), "/vsearch-1.10.2_linux_x86_64", sep="")
Sys.chmod(vsearchpath, mode = "0777", use_umask = TRUE)
}
}



if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}

if(!file.exists(file)){meep <- paste(setwd, "/", sep="")} else{meep <- ""}

cat(paste(Sys.time(), "- Clustering reads with Vsearch\n"), file=logfile, sep="", append=T)

dir.create(paste(setwd, "/Vsearch", sep=""), showWarnings=F)
filename <- gsub("(.*).fasta$", "\\1", file)

A <- system2(vsearchpath, paste(" -derep_fulllength ", meep, file, " -output ", setwd, "/Vsearch/", filename, "_drep.fasta >", setwd, "/Vsearch/temp.txt", sep=""), stdout=T, stderr=T) # dereplicate!


temp <- readLines(paste(setwd, "/Vsearch/temp.txt", sep=""))
nuc <- temp[grep("nt in", temp)]
sequ <- as.numeric(sub(".*in (.*) seqs.*", "\\1", nuc))
derep <- temp[grep("unique sequences", temp)]
derep <- as.numeric(sub("(.*)unique sequences.*", "\\1", derep))
version <- temp[grep("vsearch", temp)]
version <- sub("(.*)unique sequences.*", "\\1", version[1])

print(paste("Rading", file, collapse=""))
print(paste(sequ, "imput sequences ->", derep, "dereplicated sequences", collapse=""))

# add ";size=1"
A <- system2(vsearchpath, paste(" -derep_fulllength ", setwd, "/Vsearch/", filename, "_drep.fasta", " -output ", setwd, "/Vsearch/", filename, "_drep+1.fasta -sizeout", sep=""), stdout=F, stderr=F)



# cluster sequences!
A <- system2(vsearchpath, paste(" -cluster_fast ", setwd, "/Vsearch/", filename, "_drep+1.fasta -strand both", " -id ", id, " -msaout ", setwd, "/Vsearch/cluster_file", cmd, " >", setwd,"/Vsearch/temp.txt", sep=""), stdout=T, stderr=T) # dereplicate!


temp <- readLines(paste(setwd, "/Vsearch/temp.txt", sep=""))
clust_no <- temp[grep("Clusters: ", temp)]
clust_no <- sub("Clusters: (.*) Size.*", "\\1", clust_no)
print(paste("Clusters: ", clust_no, collapse=""))


cat(paste(version, "\n\n", sep=""), file=logfile, sep="", append=T)
cat(paste("Used fasta file: ", file, "\nNumber of imput sequences: ", sequ, "\nDereplicated: ", derep, "\nCluster: ", clust_no, "\n\n", sep=""), file=logfile, sep="", append=T)

cat(paste("VSEARCH comands:\n\n", vsearchpath, " -derep_fulllength ", meep, file, " -output ", setwd, "/Vsearch/", filename, "_drep.fasta >", setwd, "/Vsearch/temp.txt\n", 
vsearchpath, " -derep_fulllength ", setwd, "/Vsearch/", filename, "_drep.fasta", " -output ", setwd, "/Vsearch/", filename, "_drep+1.fasta -sizeout", "\n",
vsearchpath, " -cluster_fast ", setwd, "/Vsearch/", filename, "_drep+1.fasta -strand both", " -id ", id, " -msaout ", setwd, "/Vsearch/cluster_file", cmd, "  >", setwd,"/Vsearch/temp.txt", 

"\n\n", sep=""), file=logfile, sep="", append=T)


# write single fasta files from 
dir.create(paste(setwd, "/Vsearch/cluster_fasta", sep=""), showWarnings=F)

data <- readLines(paste(setwd, "/Vsearch/cluster_file", sep=""))


start <- which(data=="")
stop <- which(data==">consensus")


for (i in 1:length(start)){
sequ <- data[(start[i]+1):(stop[i]-1)]

cat(sequ, file=paste(setwd, "/Vsearch/cluster_fasta", "/", i, ".fasta", sep=""), sep="\n")
}




cat("", file=paste(setwd, "/", filename, "_cons_cluster_", threshold, ".fasta", sep=""), append=F)

# build consensus sequences!


upac <- read.csv(text=c("ID,comment,A,T,C,G,farbe
A,Adenine,1,0,0,0,F
C,Cytosine,0,0,1,0,F
G,Guanine,0,0,0,1,F
T,Thymine,0,1,0,0,F
R,A or G,0.5,0,0,0.5,T
Y,C or T,0,0.5,0.5,0,T
S,G or C,0,0,0.5,0.5,T
W,A or T,0.5,0.5,0,0,T
K,G or T,0,0.5,0,0.5,T
M,A or C,0.5,0,0.5,0,T
B,C or G or T,0,0.3,0.3,0.3,T
D,A or G or T,0.3,0.3,0,0.3,T
H,A or C or T,0.3,0.3,0.3,0,T
V,A or C or G,0.3,0,0.3,0.3,T
N,any base,0.25,0.25,0.25,0.25,T
-,gap,0,0,0,0,T"), stringsAsFactors=F)

upac[upac==0.3] <- 1/3
upac.score <- upac[,3:6] > 0

upac[upac==0.25] <- 0 # make N zero score!

letter <- c()
for (j in 1:nrow(upac.score)){
letter[j] <- paste(as.numeric(upac.score[j,]), collapse="")
}


# loop import!
d <- 6
for (d in 1:length(start)){

cat(d, file=paste(setwd, "log_cluster.txt"), sep="\n", append=T)


alignment1 <- read.fasta(paste(setwd, "/Vsearch/cluster_fasta", "/", d, ".fasta", sep=""), seqonly=T)
alignment  <- toupper(alignment1)
alignment <- unlist(strsplit(alignment, split=""))
alignment <- match(alignment, upac$ID)
alignment <- matrix(alignment, nrow=length(alignment1), ncol=nchar(alignment1[1]), byrow=T)



meep2 <- c()
for (k in 1:ncol(alignment)){
data <- alignment[,k]

colu <- upac[data, 3:6]

hey <- colSums(colu)

# shanon entropy
p <- hey / sum(hey)

if(threshold=="Majority"){
p <- p ==max(p)
} else {p <- p>=threshold}

meep2 <- c(meep2, paste(as.numeric(p), collapse=""))
}

buildsequ <- upac$ID[match(meep2, letter)]
cat(paste(">", d, "\n", paste(buildsequ, sep="", collapse=""), "\n", sep=""), file=paste(setwd, "/", filename, "_cons_cluster_", threshold, ".fasta", sep=""), append=T)


}


}




