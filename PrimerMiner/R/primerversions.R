primerversions <- function(primer, fw=T, convertInosine=T){
primer <- unlist(strsplit(primer, ""))

# nucleotide table
upac <- read.csv(text=c("ID,comment,A,T,C,G, comp
A,Adenine,1,0,0,0,T
C,Cytosine,0,0,1,0,G
G,Guanine,0,0,0,1,C
T,Thymine,0,1,0,0,A
R,A or G,0.5,0,0,0.5,Y
Y,C or T,0,0.5,0.5,0,R
S,G or C,0,0,0.5,0.5,S
W,A or T,0.5,0.5,0,0,W
K,G or T,0,0.5,0,0.5,M
M,A or C,0.5,0,0.5,0,K
B,C or G or T,0,0.3,0.3,0.3,V
D,A or G or T,0.3,0.3,0,0.3,H
H,A or C or T,0.3,0.3,0.3,0,D
V,A or C or G,0.3,0,0.3,0.3,B
N,any base,0.25,0.25,0.25,0.25,N
I,inosine,0.25,0.25,0.25,0.25,N
-,gap,0,0,0,0,-"), stringsAsFactors=F)

upac[upac==0.3] <- 1/3
upac.score <- rowSums(upac[,3:6] > 0)

if(!convertInosine){upac.score[16] <- 1}

temp <- upac.score[match(primer, upac$ID)]

val <- 1
for (i in 1:length(temp)){
val <- val*temp[i]
}

message(paste("Nuber of primer versions:", val))


# make primer versions
upac.score <- upac[,3:6] > 0
i <- 1
save <- NULL


for (i in if(fw){1:length(temp)}else{length(temp):1}){

rep_nuc <- names(upac[3:6])[upac.score[match(primer[i], upac$ID),]]

nuc <- rep(rep_nuc, val)[1:val]

if(!fw){nuc <- comp(nuc, forceToLower=F)}

if(!convertInosine & match(primer[i], upac$ID)==16){nuc <- "I"}

save <- paste(save, nuc, sep="")
save <- sort(save)
}

return(save)
}
