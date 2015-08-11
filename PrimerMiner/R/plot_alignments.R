plot_alignments <- function(files=c("List of files"), Order_names=NULL, start=NULL, end=NULL, threshold=0.1){

entropy <- function(fasta){
{
	
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
N,any base,0,0,0,0,T
-,gap,0,0,0,0,T"), stringsAsFactors=F)

alignment <- read.fasta(fasta, seqonly=T)
alignment <- strsplit(unlist(alignment), split="")
alignment <- matrix(unlist(alignment), nrow=length(alignment), ncol=length(alignment[[1]]), byrow=T)

}

meep <- c()
for (i in 1:ncol(alignment)){
data <- alignment[,i]

temp <- match(data, upac$ID)
colu <- upac[temp, 3:6]

hey <- colSums(colu)

# shanon entropy
p <- hey / sum(hey)
meep <- rbind(meep, c("ID"=i, p))

}

return(meep)
}


# Primers in one plot!!
covcol <- c("Green", "Orange", "Red", "Black")

list <- c()
for (x in 1:length(files)){
list[[x]] <- entropy(files[x])
}

plotcol <- c("Red", "Green", "Blue", "Yellow")


if (is.null(Order_names)){Order_names <- sub(".*/(.*)", "\\1", files)}
if (is.null(start)){start <- 1}
if (is.null(end)){end <- nrow(list[[1]])}


plot(-3, xlim=c(-(500/end),(abs(start-end)+2)), ylim=c(0,length(files)+5), xaxt="n", xlab="", ylab="", yaxt="n")

for (k in length(files):1){

#plot nucleotide composition
padd <- cbind(0, list[[k]][start:end,2:5])

for (i in 1:nrow(padd)){
for (m in 2:5){
rect(i-0.5, sum(padd[i,1:m-1])+3+1*k, i+0.5, sum(padd[i, 1:m])+3+1*k, col=plotcol[m-1])}
}
text(-0.5, 3.5+1*k, Order_names[k], pos=2)

}

# below berechnen
meep <- start:end

for (h in 1:nrow(padd)){
s <- 1
temp <- list[[s]][meep[h], 2:5]
if(length(files)>1){
for (s in 2:length(files)){
temp <- rbind(temp+list[[s]][meep[h], 2:5])
}
}
p <- temp / sum(temp)

#make plot

padd <- cbind(0, p*2)
for (m in 2:5){
rect(h-0.5, sum(padd[1:m-1])+1.75, h+0.5, sum(padd[1:m])+1.75, col=plotcol[m-1])}

rect(h-0.5, 0, h+0.5, 0.6, col= plotcol[which.max(p)], border=NA)
text(h, 0.3, names(as.data.frame(p))[which.max(p)], cex=0.75)

# plot nucleotide cov boxes
rect(h-0.5, 0.8, h+0.5, 1.4, col= covcol[length(which(p>threshold))], border=NA)


}

axis(side=1, at=seq(0, abs(start-end-1), 10), lab=seq(start-1, end, 10), las=2)
axis(side=1, at=seq(5, end, 10), lab=F)

}
