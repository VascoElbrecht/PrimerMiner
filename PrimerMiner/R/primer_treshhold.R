
primer_threshold <- function(fw=NULL, rw=NULL, threshold=60, file=NULL){

if(is.null(c(fw, rw))){stop("Please provide at least one table with primer evaluations! see ?evaluate_primer for details = )")}


if(!is.null(rw) & !is.null(fw)){ # compare primers
fw_data <- read.csv(fw, stringsAsFactors=F)
rw_data <- read.csv(rw, stringsAsFactors=F)

temp <- cbind("fw"=fw_data$sum, "rw"=rw_data$sum, "sum"=fw_data$sum + rw_data$sum)

if(!is.null(file)){
	write.csv(temp, file=file)
}

fail <- 0
ok <- 0
gaps <- 0

fail <- sum(temp[,3]>threshold, na.rm=T)
ok <- sum(temp[,3]<=threshold, na.rm=T)
gaps <- sum(is.na(temp[,3]))

return(unlist(list("OK"=ok, "fail"=fail, "missing"=gaps)))
} else {# evaluate single primer

if(is.null(fw)){fw <- rw}
temp <- read.csv(fw, stringsAsFactors=F)

fail <- 0
ok <- 0
gaps <- 0

fail <- sum(temp$sum>threshold, na.rm=T)
ok <- sum(temp$sum<=threshold, na.rm=T)
gaps <- sum(is.na(temp$sum))

return(unlist(list("OK"=ok, "fail"=fail, "missing"=gaps)))
}
}


#prompt(primer_threshold, "primer_treshold.Rd")
