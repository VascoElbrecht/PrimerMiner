\name{Clustering}
\alias{Clustering}
\title{
Clusters sequences by simmilarity (using Vsearch)}
\description{
You need to have VSEARCH installed for this to work!

Takes a fasta file and clusters sequences by simmilarity \code{3\%}. Before clustering the reads are  dereplicates to remove duplicates. Clustering is done with the Vsearch comand \code{cluster_fast}. Files are automatically written in the folder Vsearch. A consensus sequence of each cluster is generated within R and exported as a fasta file, whcih can be imported into Geneious.}
\usage{
Clustering(file, vsearchpath = "Vsearch", id = 0.97, cmd = "", threshold = "Majority")
}
\arguments{
  \item{file}{Imput fasta file.}
  \item{vsearchpath}{
Path to yor Vsearch program, or if it is globally executable the name of the program (usually \code{Vsearch})}
  \item{id}{
Percent simmilarity with which the reads schould be clustered. 2-3 percend are reccomended. Defult is 3 percent.}
  \item{cmd}{
Additional Vesearch commands for \code{cluster_fast} can be added here.}
  \item{threshold}{Gaps are ignored when building the consensus! The defult setting \code{Majorety} selects the most common base at each possition to build the consensus sequence. With a threshold of \code{0 - 0.25} one can define the minimum occurence of an base to be counted. A threshhold of \code{0} will lead to most ambigouties and is not recommendet.}
}
\details{
R builds its own consensus sequences from the cluster files, as Vsearch does considder gaps \code{-} when building the consensus, leading to shorter sequences.}
\value{
}
\references{
}
\author{
}
\note{
}


\seealso{
\code{\link{Merge_fasta}}
}
\examples{}