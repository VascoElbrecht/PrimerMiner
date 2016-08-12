# PrimerMiner

PrimerMiner is a R based batch sequence downloader to design and verify metabarcoding primers. Sequences for a specified marker (e.g. COI) are obtained from NCBI and BOLD and clustered into Operational taxonomic units (OTU) to reduce bias introduced by over represented sequences in the data bases.

PrimerMiner requires an internet connection for downloading sequences and was only tested on MacOSX (should work on linux as well). 

# Documentation and video tutorials

Please take a look at the PrimerMiner wiki for detailed package documentation and tutorials. We are happy to help should you run into any issues or run into problems (contact: luckylion07@googlemail.com).


## Quick guide

1) Installation
`install.packages("path_to_package_file", repos = NULL, type="source", dependencies=T)`
Load the package with `library("PrimerMiner")`. You find al the commands and an executable example in the Sample_Data!

2) Batch downloading sequences
Generate the configuration file by runing `batch_config("config.txt")`in R. Edit this file to review and change the default settings.
Create a csv table containing the groups (and their subgroups if you want to download a subset of that group) for which data should be downloaded. See "taxa_small.csv" in the folder Sample_Data.

3) Batch download sequence data
To start the batch download, run `batch_download("taxa_small.csv", "config.txt")` giving the name of the taxon table and configuration file.

4) Bild alignments
Align OTUs and extract region interesting for primer development e.g. with Geneious and MAFFT. Export the aligned sequences as a fasta file

5) Visualise alignments
With `plot_alignments(path_to_fasta_alignment_files)` you can produce plots of the alignments, to use for primer design and visualisation (see [Poster_introduction.pdf](https://github.com/VascoElbrecht/PrimerMiner/blob/master/Poster_introduction.pdf) for an example).

6) In silico Primer evaluations are now available, see the [wiki](https://github.com/VascoElbrecht/PrimerMiner/wiki/6-Primer-evaluation-(in-silico)) for more information.

## PrePrint for this package

Elbrecht & Leese (2016). https://peerj.com/preprints/2044/


## Contact

Vasco Elbrecht  - luckylion07@googlemail.com
