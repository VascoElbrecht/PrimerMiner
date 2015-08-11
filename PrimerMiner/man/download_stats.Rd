\name{download_stats}
\alias{download_stats}
\title{Summary stats for batch_download}
\description{
After using \code{batch_download} you can generate a table listing the number of sequences downloaded for each taxon from BOLD, GenBank and mitochondrial genomes.}
\usage{
download_stats(table, save = "Download_Stats_table.csv", cluster="")
}
\arguments{
  \item{table}{Supply the table used with \code{batch_download}
}
  \item{save}{Give a name under which the table should be saved as a csv.
}
  \item{cluster}{To save statistics from VSEARCH clustering, enter a docment name here. I.e. \code{"Cluster_stats.csv"}.}
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
\code{\link{batch_download}}
}
\examples{
}
