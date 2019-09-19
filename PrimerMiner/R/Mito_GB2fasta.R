Mito_GB2fasta <- function(files, marker=c("COX1", "COXi", "CO1", "COi"), add=100, rm_dup=T, pattern="\\.gb$", no_marker=T, setwd=NULL){

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}

cat(paste(Sys.time(), " - Converting Mito Genbank to fasta\n\n"), file=logfile, sep="", append=T)


cat("#mito_gb2fasta\nGBfile\tnoMito\tunique\n", file=logfile, sep="", append=T)

readin <- grep(pattern, files)



if (length(readin)==0){
	folders <- files
	readin <- c()}else{
folders <- files[-readin]
readin <- files[readin]}

no_coi_mito <- NULL

for (i in folders){
readin <- c(list.files(path=i, pattern=pattern, full.names=T), readin)
}


folder_path <- paste(folders, "/", sep="")


for (k in 1:length(readin)){

# extract fasta file!
gb <- readLines(readin[k])
sequence_exp <- NULL

range <- grep("//", gb)
range <- c(1, range)
for (i in 1:length(range[-1])){
temp <- gb[range[i]:range[i+1]]

query <- paste("/(gene|CDS|product)=\"(", paste(marker, sep="", collapse="|") ,")\"", sep="")
cds_pos <- grep("( gene   | CDS   | rRNA   )", temp)

COX1 <- cds_pos[grep(query, temp[cds_pos+1], ignore.case=T)[1]]
if(is.na(COX1)){COX1 <- cds_pos[grep(query, temp[cds_pos+2], ignore.case=T)[1]]} #try second line!
COX1 <- temp[COX1]
if(grepl("join", COX1)){COX1 <- c(NA, NA)} else{
COX1 <- gsub("[[:alpha:]<> ()]", "", COX1)
COX1 <- as.numeric(c(gsub("(.*)\\.\\..*$", "\\1", COX1), gsub(".*\\.\\.(.*)$", "\\1", COX1)))}

if (!is.na(COX1[1])&!is.na(COX1[2])&length(grep("ORIGIN      ", temp))!=0){
sequ <- temp[(grep("ORIGIN      ", temp)+2):length(temp)-1]
sequ <- paste(gsub("[ 1234567890]", "", sequ), collapse="")

sequ <- substr(sequ, COX1[1]-add, COX1[2]+add)
sequ_name <- gsub(".* (.*)$", "\\1", temp[grep("ACCESSION ", temp)])


sequence_exp <- c(sequence_exp, paste(">", sequ_name, "\n", sep=""), paste(sequ, "\n", sep=""))} else { # export mito with no marker
no_coi_mito <- c(no_coi_mito, gsub(".* (.*)$", "\\1", temp[grep("ACCESSION ", temp)]))
if (no_marker){dir.create("mito_no_marker", showWarnings=F)
Fname <- sub(pattern, "\\1", readin[k])
Fname <- sub(".*/(.*)", "\\1", Fname)
#torder <- readLines(logfile)[grep("#mito_data: Folder .*", readLines(logfile))][1]
#torder <- sub("#mito_data: Folder (.*)", "\\1", torder)
Fname1 <- paste("mito_no_marker/", setwd, if(!is.null(setwd)){"_"}, Fname, ".gb", sep="")

if(temp[1]=="//") ( temp <- temp[-c(1,2)])

cat(temp, file=Fname1, sep="\n", append=T)}
}
}


if (!is.null(sequence_exp)){
if (rm_dup){
dup <- !duplicated(sequence_exp[seq(2, length(sequence_exp), 2)])
sequence_exp <- sequence_exp[rep(dup, each = 2)]}

cat(sequence_exp, file=paste(sub(pattern, "\\1", readin[k]), ".fasta", sep=""), sep="") # overwrite old files

# check number of unique sequences
uni_mito_coi <- length(unique(sequence_exp[seq(2, length(sequence_exp), 2)])) } else {uni_mito_coi <- 0}


message(paste("Converted ", length(range[-1])," mitogenome", if(length(range[-1])==1){""}else{"s"}," of ", readin[k], " to ", uni_mito_coi, " unique ", marker[1], " sequence", if(uni_mito_coi==1){""}else{"s"},".", sep=""))
cat(paste(readin[k],"\t", length(range[-1]), "\t", uni_mito_coi, "\n", sep=""), file=logfile, sep="", append=T)
}

cat("#mito_gb2fasta_end\n\n", file=logfile, sep="", append=T)

if (!is.null(no_coi_mito)){
cat(paste("Mitogenomes in which no marker sequences were detected:", paste(no_coi_mito, collapse=", "), "\n\n"), file=logfile, sep="", append=T)
message(" ")
message(paste("Mitogenomes in which no marker sequences were detected:", paste(no_coi_mito, collapse=", ")))}
}
