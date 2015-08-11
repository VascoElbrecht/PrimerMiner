\name{plot_alignments}
\alias{plot_alignments}
\title{
Plot all provided alignments into one overview to design primers.
}
\description{
Plot the base proportions for each alignment at each position to manually find the best primers. All alignment files have to be in the fasta format and should be already globally aligned (same length).}
\usage{
plot_alignments(files = c("List of files"), Order_names = NULL, start = NULL, end = NULL, threshold = 0.1)
}
\arguments{
  \item{files}{
lift of sastafiles to plot (include full file path!).}
  \item{Order_names}{
If you want to include your own names for each fasta file, write them here. Otherwise the fasta filename will be used.}
  \item{start}{
Select from which base pair position start plotting (default = 1)}
  \item{end}{
Select from which base pair position end plotting plotting (default = plot complete alignment)}
  \item{threshold}{
Minimum \% abundance a base needs to have to be counted in the summary plot (default = 0.1)}
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
}
\examples{
}
