\name{primer_treshold}
\alias{primer_treshold}
\title{
Analyse evaluate_primer() tables
}
\description{
Takes tables generated with \code{\link{evaluate_primer}} and sums up the number of primers which worked and failed based on a defined \code{threshold}. When two files are given, the mismatch penalties from both primers are added together before sequences are filtered based on the threshold (scores can be saved with \code{file}). It is recommended to use a higher threshold when evaluation primer pairs than when looking at single primers.
}
\usage{
primer_treshold(fw = NULL, rw = NULL, threshold = 60, file = NULL)
}
\arguments{
  \item{fw}{
Give the path to the first primer table (e.g. forward).}
  \item{rw}{
Give the path to the second primer table (e.g. reverse).}
  \item{threshold}{
All sequences with a penalty score of above the \code{threshold} will be considered "not working". Sequences with a score of NA (because of missing data, gaps or N in the sequence) are counted as "missing".}
  \item{file}{
When processing tables from two primers, the summed up penalty scores can be written in a file (supply file name in this variable).}
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
\code{\link{evaluate_primer}}
}
\examples{
}