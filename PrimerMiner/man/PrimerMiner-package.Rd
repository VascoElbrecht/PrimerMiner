\name{PrimerMiner-package}
\alias{PrimerMiner-package}
\alias{PrimerMiner}
\docType{package}
\title{Generate and visualise degenerated Primers from DNA sequences}
\description{This packges can download DNA sequences of the targed marker from the NCBI and BOLD database, to design degenerate primers for specific taxonomic groups. For Sequence clustering Vsearch is neded, and for sequence assembly Geneious is required.}
\details{
\tabular{ll}{
Package: \tab PrimerMiner\cr
Type: \tab Package\cr
Version: \tab v0.01\cr
Date: \tab 2015-04-06\cr
License: \tab GPL-3\cr
}
PrimerMiner can download and cluster sequence data for specific taxonomic groups, which can then be used to generate optimesed barcoding primers.

Target taxons are supplied in a csv file using scientific names (usaually on Order and Familiy level). In the first coloum, the Target taxa wich should be downloaded and clustered are supplied. For each Taxon, sequence processing is done seperately. In the Second coloum, the sequence search can be narriwed down to certain taxa (for example Families) within the taxonomic group. See sample files for mor information.

For downloading, clipping (optional) and clustering, diferent functions are available. However, when taxa are supplied in a taxon file, sequences can be downloaded and processed in batch. Use batch_config to create a configuration txt file, in wich individual parameters can be adjusted, and batch_download to batch download and process the files.
Besides statistic and log files, R will provide you with the clustered reads in the fasta file format, which can be directly imported into Geneious for assembly and primer design.

Base composition of target sites can be visualised using the XXXX function.

Please see the publication and package PDF manual for further information. The latest release of PrimerMiner is available on GitHub. Please report bugs and requests also here. 
}
\author{Vasco Elbrecht
Maintainer: Vasco Elbrecht <luckylion07@googlemail.com>
}
\references{GeneStream.de}
\keyword{ package }