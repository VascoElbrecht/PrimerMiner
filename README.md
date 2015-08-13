# PrimerMiner
PrimerMiner is a R based package, which allows you to batch download a specific barcoding gene (like COI, 16S) from NCBI and BOLD databases. The sequences are then clustered into Operational taxonomic units (OTU) using Vsearch and the consensus sequence of each cluster saved in a fasta file. The OTUs in wich each species is represented by only ~one sequence, can then be used in other programs for universal primer development, or aligned and visualised with PrimerMiner for finding good primers.

Please take a look at the [Poster_introduction.pdf](https://github.com/VascoElbrecht/PrimerMiner/blob/master/Poster_introduction.pdf) to learn more about the ideas behind PrimerMiner as well as a quick start tutorial. You can download the poster by clicking on "RAW" on the top right.

**BETA RELEASE**: This is a very early beta release! The software works, but I'm still optimising and adding documentation!

**Only tested on MacOSX** In theory this works on linux, but PrimerMiner is currently only developed and tested on MacOSX. 

## Known issues
- Not tested on Linux!
- On Mac, don't run Primer Miner from the download folder (causes a problem to generate the sumary statistiks)
- Lacks good documentation, I'm working on it
- Found another problem? Please send me your configuration file, taxa table and which operating system you are using and I will take a look (luckylion07@googlemail.com)!

## Installing PrimerMiner

PrimerMiner is package for R, wich can be downloaded on [r-project.org](https://www.r-project.org/)

You can download the latest release of [PrimerMiner](https://github.com/VascoElbrecht/PrimerMiner/releases) and install it with the following commands in R:

`install.packages(c("bold", "XML", "rentrez", "seqinr", "devtools"), dependencies=T)`
`install.packages("path_to_file", repos = NULL, type="source")`

Load the package with `library("PrimerMiner")`. You find al the commands and an executable example in the Sample_Data!

## Batch downloading sequences

Setting for marker, downloading and clustering can be made in the configuration file, wich can be created by running `batch_config("config.txt")`in R. The txt file can then be edited and settings adjusted.

By default data is downloaded for the standard COI barcoding marker. If you would like to download sequences for another marker change `Marker = c("16S", "large subunit ribosomal RNA", "16S ribosomal RNA", "l-rRNA")'. If you are downloading data for a barcoding marker which is not present on BOLD, make sure to set ‘download_bold = T‘.


In a second file you have provide information about the taxa for which sequences should be downloaded. This is done in a csv table file like in the Sample_Data. You can specify taxonomic groups with latin names, and also download a subset within a group (like certain families within the order Coleoptera). Please make sure to keep the taxonomic scope broad, on order or family level. On a genus or species level PrimerMiner does not work well, due to lack of sequences in databases and multiple hits in BOLD for the same query.

To start the batch download, run `batch_download("taxa_small.csv", "config.txt")` giving the name of the taxon table and configuration file.

PrimerMiner might take several hours to download and process the data for big datasets! Please don't touch anything and let R do it's thing.

## What to do with the PrimerMiner output?

PrimerMiner Uses [Vsearch](https://github.com/torognes/vsearch) to cluster sequences and generate consensus sequences of the OTUs. 





## References



## Contact






