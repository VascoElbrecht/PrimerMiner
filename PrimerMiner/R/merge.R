# merge fasta

Merge_fasta <- function(files, save, clip_left=26, clip_right=26, pattern="\\.fasta$", minlength=50, overwrite=T, setwd=NULL) {

if(!is.null(setwd)){logfile <- paste(setwd, "/log.txt", sep="")} else {logfile <- "log.txt"}


cat(paste(Sys.time(), "- Merging fasta files\n\n"), file=logfile, sep="", append=T)


readin <- grep(pattern, files)

if (length(readin)==0){
	folders <- files
	readin <- c()}else{
folders <- files[-readin]
readin <- files[readin]}

cat(paste("Reading in files matching ", pattern, ":", "\nFolders: ", paste(folders, collapse=", "), "\nFiles: " , paste(readin, collapse=", "), "\n", "Clipping: Left ", clip_left, " bp, Rigth ", clip_right, " bp\n\n", sep=""), file=logfile, sep="", append=T)
message(paste("Clipping: Left ", clip_left, " bp, Rigth ", clip_right, " bp", sep="", collapse=""))

if (overwrite){
if(file.exists(save)){unlink(save)
message(paste("deleted old fasta file: ", save, colapse="", sep="" ))
cat(paste("deleted old fasta file: ", save, "\n\n", colapse="", sep="" ), file=logfile, sep="", append=T)}
}

for (i in folders){
readin <- c(list.files(path=i, pattern=pattern, full.names=T), readin)
}
message("Files writen in fasta")
message(readin)
message(" ")

cat(paste("Matching files which were written into ", save, ":", "\n", paste(readin, collapse=", ", sep=""), "\n\n"), file=logfile, sep="", append=T)

for (k in 1:length(readin)){ # read in indiv files
fasta <- read.fasta(readin[k], as.string=T)

empty <- c()

for (j in 1:length(fasta)){ # apply clipping
fasta[j][[1]][1] <- substr(fasta[j][[1]][1], clip_left+1, nchar(fasta[j][[1]][1])-clip_right)
if (nchar(fasta[j][[1]][1])<minlength) {empty <- c(empty, j)}
}
if (is.null(empty)){write.fasta(fasta, names=names(fasta), file.out=save, open="a")} else {write.fasta(fasta[-empty], names=names(fasta)[-empty], file.out=save, open="a")}
}
if (!is.null(empty)){cat(paste(length(empty), " sequences discarded as they have a length below ", minlength, " bp", "\n\n", sep=""), file=logfile, sep="", append=T)}
}


