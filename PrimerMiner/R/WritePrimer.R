WritePrimer <- function(file, sequence, name="NA", type="Ref"){

# chack if file exists
if(!file.exists(file)){
cat("Type\tprimer\tN\tcov\tsd\theadline\tfile_used", file=file, sep="\n")
print(paste("File ", file, " does not exist and was therfore created!", sep=""))
}

if (type=="Line"){
out <- paste("Line", "NA", "NA", "NA", "NA", "NA", "NA", sep="\t")
cat(out, file=file, append=T, sep="\n")
print(paste("A line indicator was writen to the file ", file, "!", sep=""))
}else{
sequence <- toupper(sequence)
out <- paste("Ref", sequence, "NA", "NA", "NA", name, "NA", sep="\t")
cat(out, file=file, append=T, sep="\n")# write file
print(paste("The primer  ", sequence, " was added to the file ", file, "!", sep=""))}
}
