\name{selectivetrim}
\alias{selectivetrim}
\title{Trimm Primer sequences and remove gaps}
\description{
R
}
\usage{selectivetrim(read, write, trimL = 25, trimR = 26, gaps = 0.10, minsequL = 100)}
\arguments{
  \item{read}{Path to imput fasta file.}
  \item{write}{Save processed fasta file as "file.name".}
  \item{trimL}{Apply left sited selective trimming.}
  \item{trimR}{Apply right sited selective trimming.}
  \item{gaps}{Remove all gaps with not at least X\% sequences.}
  \item{minsequL}{Sequences shorter than minsequL are discarded.}
}
\details{
Generated alignments still contain gaps and potentially still primer sequences if no primer trimming was applied when \code{\link{batch_download}} sequences.

With \code{gaps=0.10} most gaps can be removed. Doble check that the sequence length of the alignemnt is as expected and adjust the treshhold accordingly. 

Additionally primer trimming from the left and right can be applied, which does not affect the actual amplified region. The length of the primer has to be included as extra. R applies X bp clipping in X*2 in base pairs from the left or right. 
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
\code{\link{batch_download}}, \code{\link{plot_alignments}}
}
\examples{
}