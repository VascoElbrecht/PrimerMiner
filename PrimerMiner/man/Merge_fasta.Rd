\name{Merge_fasta}
\alias{Merge_fasta}
\title{
Merges several fasta files into one (with clipping)
}
\description{
Read in single fasta files or Folders containing fasta files and merge them into one single fasta file. With defult settings, 26 bp are clipped from the right and left, to remove potential primer seuences in the data (\code{clip_left} and \code{clip_right}). Sequences of \code{minlength} 50 bp are not writen in the fasta file. 
}
\usage{
Merge_fasta(files, save, clip_left = 26, clip_right = 26, pattern = "\\.fasta$", minlength = 50, overwrite = T)
}
\arguments{
  \item{files}{
A vector containing fasta file names of folder names which contain fasta files. Mixing fasta files and folders is possible. To ajust files which are being detected take a look at \code{pattern}. 
}
  \item{save}{
Name of the document in which the files should be saved. For example \code{All_plecoptera.fasta}}
  \item{clip_left}{
Number of bp clipped from the left. Set to \code{0} if you don't wich to apply clipping to the sequences.}
  \item{clip_right}{
Number of bp clipped from the right. Set to \code{0} if you don't wich to apply clipping to the sequences.}
  \item{pattern}{
Pattern defining which fasta files are detected when folders are given. By defult all \code{.fasta} files are detected with \code{"\\.fasta$"}. To only load files downloaded from one source add \code{"BOLD\\\\.fasta$"}, \code{"GB\\\.fasta$"} or \code{"mito\\\.fasta$"}}
  \item{minlength}{
Define minimum length sequences should have to be written in the output fasta file. Minimum length is 1 bp.}
  \item{overwrite}{
Defult is TRUE, which overrides the file if \code{save} already exists. If \code{overwrite=TRUE} reads are appended if the file exists. This can be useful to merge sequences from BOLD and GenBank with clipping apllied, and than add mitochondrial reads without clipping. 
}
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
\code{\link{Download_BOLD}}, \code{\link{Download_GB}}, \code{\link{Download_mito}} and \code{\link{Clustering}}
}
\examples{}