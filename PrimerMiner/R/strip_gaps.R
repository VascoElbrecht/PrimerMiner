# 160905 - strip gaps

selectivetrim <- function(read, write, trimL=25, trimR=26, gaps=0.10){

data <- read.fasta(read, as.string=T, forceDNAtolower=F) # read file
data <- unlist(data)
sequnames <- names(data)


left <- as.vector(nchar(sub("(-*).*", "\\1", data)))
right <- as.vector(nchar(sub(".*[ACGTRYSWKMBDHVN](-*)", "\\1", data)))
right <- as.vector(nchar(data[1])) - right
lengtbefore <- as.vector(nchar(data[1])) 


data <- strsplit(data, "") # split
data <- matrix(unlist(data), length(data), byrow=T) # make matrix


if(is.numeric(gaps)){ # remove gaps
exp <- NULL

for (i in 1:ncol(data)) {
temp <- !left>=i # left
temp[which(temp)] <- right[which(temp)] >= i # check right
temp <- data[,i][temp]

exp <- rbind(exp, cbind(length(temp), sum(temp!="-")))
}

keep <- exp[,2]/exp[,1]
data <- data[,(keep>=0.10)]

message(paste("Number of bases", lengtbefore, "before, removing gaps with >=", gaps, "gaps -> alignment length now", ncol(data), "bp", collapse=""))
}

# apply selective trimming!

if (trimL>0){
meep <- NULL
for (i in 1:nrow(data)){
meep[i] <- paste(data[i,], collapse="")
}

left <- as.vector(nchar(sub("(-*).*", "\\1", meep)))
trim <- trimL*2>left # trim only up to primer bind
trim <- which(trim)


trimto <- trimL +left[trim] 
trimto[trimto>trimL*2] <- trimL*2

#cbind(left[trim], trimto, trim)

for (i in 1:length(trim)){
data[trim[i], left[trim[i]]:trimto[i]] <- "-"
}

data <- data[, -(0:trimL)]

message(paste("Left sided selective trimming applied with", trimL, "bp", collapse=""))
}




# selcetive trimming on the rigth!

if (trimR>0){
meep <- NULL
for (i in 1:nrow(data)){
meep[i] <- paste(data[i,], collapse="")
}

sequlength <- nchar(meep[1])

rigth <- as.vector(nchar(sub(".*[ACGTRYSWKMBDHVN](-*)", "\\1", meep)))
trim <- trimR*2>rigth # trim only up to primer bind
trim <- which(trim)


trimto <- trimR + rigth[trim] 
trimto[trimto>trimR*2] <- trimR*2


trimto <- trimto-1
#rigth <- rigth-1

#cbind(left[trim], trimto, trim)

for (i in 1:length(trim)){
data[trim[i], (sequlength-trimto[i]):(sequlength-rigth[trim[i]])] <- "-"
}

data <- data[, -((sequlength+1-trimR): sequlength)]

message(paste("Right sided selective trimming applied with", trimR, "bp", collapse=""))
}


# prepare sequence for saving
meep <- NULL
for (i in 1:nrow(data)){
meep[i] <- paste(data[i,], collapse="")
}

write.fasta(sequences=as.list(meep), names=as.list(sequnames), write)

}



