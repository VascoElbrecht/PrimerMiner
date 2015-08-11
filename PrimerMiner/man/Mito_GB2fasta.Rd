\name{Mito_GB2fasta}
\alias{Mito_GB2fasta}
\title{
Extract marker from GenBank file and save as fasta
}
\description{
With \code{\link{Download_mito}} we ahve downloadaed mitochondrial genomes for our target taxa. Now, we need to extract the bacoding gene of interest and save the genes as a fasta file. With the defult settings, the COi marker is searched and 100 bp are added to the extracted COi gene, to account for annotation margins.
}
\usage{
Mito_GB2fasta(files, marker = c("COX1", "CO1", "COi", "COiX"), add = 100, rm_dup = T, pattern = "\\.gb$")
}
\arguments{
  \item{files}{
Path to individual GenBank files or folder containing GenBank files ending with \code{.gb}.}
  \item{marker}{
Name of the marker in the GenBank file. Capitalisation does not madder. If your target Gene has diferent common spellings include them all like: \code{c("MyGene", "MiGene", "MGene")}}
  \item{add}{
If possible add \code{add=100} bp on the left and right on the detected marker sequence, to account for lengthvariation in the annotations.}
  \item{rm_dup}{
When \code{TRUE} the marker sequences are derepilcated. Many mitochondrial Genomes are present in GenBank twice. As we will be using the mitochondrial genomes as a reference to mapp all our reads, It's highly recommended to keep this option on!}
  \item{pattern}{
Specify which .gb files are detected. At defult all documends endeing with .gb in a folder are read.}
  \item{no_marker}{
if TRUE, write mitochodrial genomes which contain no hits for \code{marker} as a GenBank file for inspection in Geneious in the folder \code{"mito_no_marker"}.}
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