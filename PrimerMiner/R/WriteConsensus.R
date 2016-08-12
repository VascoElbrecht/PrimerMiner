WriteConsensus <- function(file, export=file, threshold=0){

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


data <- read.csv(file, sep="\t", stringsAsFactors=F)
data <- data[data$Type=="primer",]


primer <- data$primer
primer <- strsplit(primer, split="")
pleng <- length(primer[[1]])
pnumber <- length(primer)


consensus <- c()
export <- upac[0, 3:6]

for (m in 1:pleng){
temp <- c(0,0,0,0)
for (n in 1:pnumber){
	meep <- which(primer[[n]][m]==upac$ID)
	temp <- c(temp+upac[meep, 3:6])
}
export <- rbind(export, temp)
}

temp2 <- upac[,3:6]>0
compare <- c()
for (j in 1:nrow(temp2)){
compare[j] <-  paste(temp2[j,], collapse="")
}

export <- export*100/pnumber
export[export == 0] <- -1
export <- export>=threshold


letter <- c()
for (j in 1:nrow(export)){
letter[j] <- which(paste(export[j,], collapse="")==compare)
}


# check if file exists
if(!file.exists(file)){
cat("Type\tprimer\tN\tcov\tsd\theadline\tfile_used", file=file, sep="\n")
message(paste("File ", file, " does not exist and was therfore created!", sep=""))
}

sequence <- paste(upac$ID[letter], sep="", collapse="")
out <- paste("Ref", sequence, "NA", "NA", "NA", "Consensus", "NA", sep="\t")
cat(out, file=file, append=T, sep="\n")# write file
message(paste("The consensus sequence  ", sequence, " was added to the file ", file, "!", sep=""))

}

