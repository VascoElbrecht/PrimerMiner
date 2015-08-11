\name{Download_BOLD}
\alias{Download_BOLD}
\title{Downloas sequences from BOLD}
\description{Will download sequences from the BOLD database and save them as fasta. Can download data for multiple taxa if \code{taxon} contains several taxa names. Will save files in a folder if given in \code{folder}. Statistics are written in the log.txt file.}
\usage{
Download_BOLD(taxon, folder = NULL)
}
\arguments{
  \item{taxon}{Contains taxon same which sequence data should be obtained for using the \code{BOLD} API. Can download sequences for several taxa if several names are given in a vecotor \code{c("Taxon1", "Taxon2")}.}
  \item{folder}{If a folder name is given, it will be created and downloaded files saved as a fasts file in the folder. If no folder is given, data is saved in the current directory.}
}
\details{Please make sure that the \code{taxon} can actually be found in BOLD and has an unique name. Might not work proppaly on a Genus or Species level. We recommend to search only on a Family or Order level!}
\value{
}
\references{}
\author{}
\note{}
\seealso{\code{\link{Download_GB}}, \code{\link{Download_mito}}}
\examples{
\dontrun{
# To download all sequences of the order Plecoptera run (might take a few minutes)
Download_BOLD("Plecoptera")

# Download data for selected families of the order Hirudinea into the folder "Hirudinea"
Download_BOLD(c("Erpobdellidae", "Glossiphoniidae", "Haemopidae", "Hirudinidae", "Piscicolidae", "Salifidae"), folder="Hirudinea")
}
}