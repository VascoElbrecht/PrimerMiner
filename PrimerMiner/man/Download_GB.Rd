\name{Download_GB}
\alias{Download_GB}
\title{Downloas sequences from GenBank nuccore}
\description{Will download sequences from the GenBank Nucleotide database and save them as fasta. Can download data for multiple taxa if \code{taxon} contains several taxa names. Will save files in a folder if given in \code{folder}. Statistics are written in the log.txt file.}
\usage{
Download_GB(taxon, folder = NULL, marker = c("COi", "CO1", "COXi", "COX1"), maxlength = 2000, custom_query = NULL)
}
\arguments{
  \item{taxon}{
Contains taxon same which sequence data should be obtained for using the \code{rentrez} R package. Can download sequences for several taxa if several names are given in a vecotor \code{c("Taxon1", "Taxon2")}.}
  \item{folder}{
If a folder name is given, it will be created and downloaded files saved as a fasts file in the folder. If no folder is given, data is saved in the current directory.}
  \item{marker}{
Provide commonly used abriviations for the barcdoing marker you are looking for. Defult is set up for the Coytocrome oxidase subunit 1: \code{c("COi", "CO1", "COXi", "COX1")}.}
  \item{maxlength}{
Set the maximal sequence length of the found sequences. In the defult settings this is limited to 2000 bp sequence length, as the COi gene is only ~1600 bp long}
  \item{custom_query}{
If you like to use a custom query NCBI, you can insert it here. If a custom search query is used, values from \code{marker} and \code{maxlength} are ignored and the \code{taxon} is placed right in fromnt of your \code{custom_query}}
}
\details{
}
\value{
}
\references{
}
\author{
}
\note{
}

\seealso{
\code{\link{Download_BOLD}}, \code{\link{Download_mito}}
}
\examples{
\dontrun{
# To download all sequences of the order Plecoptera run (might take a few minutes)
Download_GB("Plecoptera")

# Download data for selected families of the order Hirudinea into the folder "Hirudinea"
Download_GB(c("Erpobdellidae", "Glossiphoniidae", "Haemopidae", "Hirudinidae", "Piscicolidae", "Salifidae"), folder="Hirudinea")
}
}