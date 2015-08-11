\name{Download_mito}
\alias{Download_mito}
\title{Downloade mitochondrial genomes in the GenBank format}
\description{
Will download mitochondrial genomes from the GenBank Nucleotide database and save them as GenBank files. Can download data for multiple taxa if \code{taxon} contains several taxa names. Will save files in a folder if given in \code{folder}. Statistics are written in the log.txt file.

GenBank files can be converted to fasta with XXXX.}
\usage{
Download_mito(taxon, folder = NULL, minlength = 2001, maxlength = 80000, custom_query = NULL)
}
\arguments{
  \item{taxon}{
Contains taxon same which sequence data should be obtained for using the \code{rentrez} R package. Can download sequences for several taxa if several names are given in a vecotor \code{c("Taxon1", "Taxon2")}.}
  \item{folder}{
If a folder name is given, it will be created and downloaded files saved as a fasts file in the folder. If no folder is given, data is saved in the current directory.}
  \item{minlength}{
Set the maximal sequence length of the found sequences. Is set to 2001 bp to avaoid overlap with the settings in \code{Download_GB}}
  \item{maxlength}{
Set the maximal sequence length of the found sequences. In the defult settings this is limited to 80000 bp sequence length, as the mitochondrial genomes of animals are usually not longer. Adjust if needed!}
  \item{custom_query}{
If you like to use a custom query NCBI, you can insert it here. If a custom search query is used, values from \code{marker} and \code{maxlength} are ignored and the \code{taxon} is placed right in fromnt of your \code{custom_query}}
}
\details{
}
\value{}
\references{
}
\author{
}
\note{
}
\seealso{
}
\examples{
}
