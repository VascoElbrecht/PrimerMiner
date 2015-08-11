# PrimerMiner
PrimerMiner is a R based package, which allows you to batch download a specific barcoding gene (like COI, 16S) from NCBI and BOLD databases. The sequences are then clustered into Operational taxonomic units (OTU) using Vsearch and the consensus sequence of each cluster saved in a fasta file. The OTUs in wich each species is represented by only ~one sequence, can then be used in other programs for universal primer development, or aligned and visualised with PrimerMiner for finding good primers.

Please take a look at the [Poster_introduction.pdf](https://github.com/VascoElbrecht/PrimerMiner/blob/master/Poster_introduction.pdf) to learn more about the ideas behind PrimerMiner as well as a quick start tutorial. You can download the poster by clicking on "RAW" on the top right.

**BETA RELEASE**: This is a very early beta release! The software works, but I'm still optimising and adding documentation!

## Installing PrimerMiner

PrimerMiner is package for R, wich can be downloaded on [r-project.org](https://www.r-project.org/)



You can download the latest release of PrimerMiner by running the following commands in R

´´´
install.packages(c("bold", "XML", "rentrez", "seqinr", "devtools"), dependencies=T)
library("devtools")
install_github("VascoElbrecht/PrimerMiner/PrimerMiner")´´´


Alternatively you can download the latest release of PrimerMiner (which includes the sample files) and install the package with

install.packages("path_to_file", repos = NULL, type="source")


## References





