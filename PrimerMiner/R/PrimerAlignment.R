PrimerAlignment <- function(file, headline="", mean=T, spaceleft=-9){

data <- read.csv(file, sep="\t", stringsAsFactors=F)

lines <- which(data$Type=="Line")
dim <-nrow(data)-length(lines)*0.5
pleng <- length(strsplit(data[-lines,][1,2], split="")[[1]])

letters <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "I")


# make plot
plot(-2,-6, xlim=c(spaceleft, if(mean){pleng+4}else{pleng}), ylim=c(1, dim), xaxt="n", yaxt="n", xlab="", ylab="")

for (i in 1:nrow(data)){
temp <- data[i,]

if (temp$Type=="Line"){
	lines(c(0, pleng)+0.5, c(dim, dim)+ 0.25)
	dim <- dim-0.5
} else {
text(0, dim, temp$headline, adj=1)
primersequ <- unlist(strsplit(temp$primer, split=""))
for (k in which(primersequ%in%c(letters))){
	rect(k-0.5, dim-0.5, k+0.5, dim+0.5, col="Lightgray", border=NA)}
text(1:pleng, dim, primersequ)
if(mean){if(!is.na(temp$cov)){text(pleng+5, dim, temp$cov, adj=1)}}

dim <- dim-1
}
}
mtext(headline, 3, font=2, adj=0)

}
