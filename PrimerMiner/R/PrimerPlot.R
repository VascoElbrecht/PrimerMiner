PrimerPlot <- function(file, threshold=10, plot=T, otulim="auto", headline=file, export="off", rmN=T, covtab="off", pch=21){


# load upac table
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


if(!is.data.frame(file)){ # skip if not fasta file

data <- readLines(file)
data <- data[seq(2, length(data), 2)]

plength <- length(unlist(strsplit(data[1], "")))

# make matrix with single bases
temp <- strsplit(data, "")
data <- matrix(unlist(temp), nrow=length(data), ncol=plength, byrow=T)

# calculate over all coverage!
temp <- data

temp[temp=="-"] <- 0 # "-" is not calculated in the coverage
if(rmN){temp[temp=="N"] <- 0} else {temp[temp=="N"] <- 1}
temp[temp!=0] <- 1

temp <- matrix(as.numeric(temp), nrow=nrow(data), ncol=plength, byrow=F)

coverage <- colSums(temp)


# calculate indiv proportion of bases!
# count N as 0!

upacScore <- upac
if(rmN){upacScore[15,3:6] <- 0} # make N zero poins!

Texport <- upacScore[0, 3:6]
for (i in 1:ncol(data)){

colu <- upac[0, 3:6]
for(k in 1:nrow(data)){
	temp <- which(data[k,i]==upacScore$ID)
	colu <- rbind(colu, upacScore[temp, 3:6])
	
}
Texport <- rbind(Texport, colSums(colu))
}

names(Texport) <- names(upacScore[0, 3:6])} else { # if data frame was supplied
	Texport <- file
	coverage <- rowSums(Texport)
}


#####
percent <- Texport*100/coverage
rowSums(percent)

if(otulim=="percent"){otulim <- 100
coverage <- coverage*100/max(coverage)	
}

# write proportion table
if (covtab!="off"){
write.csv(Texport, file=paste(covtab, "_", gsub("/", "", headline), ".csv", sep=""))
}


# calculate primer sequence!!!

ptext <- percent
ptext[ptext == 0] <- -1
ptext <- ptext>=threshold # consensus threshold!

temp2 <- upac[,3:6]>0

compare <- c()
for (j in 1:nrow(temp2)){
compare[j] <-  paste(temp2[j,], collapse="")
}

letter <- c()
for (j in 1:nrow(ptext)){
letter[j] <- which(paste(ptext[j,], collapse="")==compare)
}

########
# write primer data in file!
if(export!="off"){

if(!file.exists(export)){
cat("Type\tprimer\tN\tcov\tsd\theadline\tfile_used", file=export, sep="\n")}



out <- paste("primer", paste(upac$ID[letter], collapse=""), nrow(data), round(mean(coverage), 2), round(sd(coverage), 2), headline, file, sep="\t")

cat(out, file=export, append=T, sep="\n")# write file


}

#############
# make plot!

if(plot){

if(otulim=="auto"){otulim <- max(coverage)*1.1}

p <- percent

par(mar=c(1,4.7,2,0.7)) # down, left, top, rigth
plot(-2,0, ylim=c(-30,220), xlim=c(1,nrow(p)), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")

for (i in 1:nrow(p)){
rect(i-0.5, 0, i+0.5, p[i,1], col="Red") #A
rect(i-0.5, p[i,1], i+0.5, rowSums(p[i,1:2]), col="Green")
rect(i-0.5, rowSums(p[i,1:2]), i+0.5, rowSums(p[i,1:3]), col="Blue")
rect(i-0.5, rowSums(p[i,1:3]), i+0.5, rowSums(p[i,1:4]), col="Yellow")
if (upac$farbe[letter][i]){
	rect(i-0.5, -10, i+0.5, -30, col="Lightgray", border=NA)}
}

rect(0.5, 120, i+0.5, 220)
rect(0.5, 0, i+0.5, 100)

points(120+coverage*100/otulim, pch= pch)
lines(120+coverage*100/otulim)

text(1:nrow(p), -20, upac$ID[letter])

axis(2, at=c(0,20,40,60,80,100), las=1, cex.axis=0.7) # base abundance
axis(2, at=c(0,20,40,60,80,100)+120, las=1, labels=round(c(0, otulim*0.2, otulim*0.4, otulim*0.6, otulim*0.8, otulim), digits=0), cex.axis=0.7) # base abundance

text(0.25, 206, paste("N = ", nrow(data), ", Mean = ", round(mean(coverage), 2), ", SD = ", round(sd(coverage), 2), sep=""), pos=4) # add coverage text

mtext(headline, side=3, at=c(0.5), adj=0, font=2)
mtext("Number of OTUs", side=2, at=c(170), las=0, line=2.2, cex=0.7)
mtext("Nucleoptide\nproportion [%]", side=2, at=c(50), las=0, line=2.2, cex=0.7)
}
# plot end 

}