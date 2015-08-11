\name{PrimerPlot}
\alias{PrimerPlot}
\title{Plot primer composition, coverage and sequence}
\description{
Reads a fasta file containing several sequences from the primer binding region (Ideally OTUs). Will create a plot showing the base composition for each position in the primer binding side, as well as the over all coverage and degenerated primer sequence.
}
\usage{PrimerPlot(file, threshold = 10, plot = T, otulim = "auto", headline = file, export = "off", rmN = T)}
\arguments{
  \item{file}{Given fasta file with multiple sequences of the primer binding region. Make sure the sequences are aligned and gaps are filled with "-".}
  \item{threshold}{Threshold at wich bases are considered for the degenerated primer sequence. The default of 10 percent is usually fine.}
  \item{plot}{Indicate if a plot should be generated. Default i TRUE}
  \item{otulim}{Set the dimensions of the coverage plot. Can be a number, which is useful if plots for forward and reverse primer should have the same dimensions. "Percent" will show the axis in percent, \code{"Auto"} will scale the plot automatically using the coverage data.}
  \item{headline}{Specify a headline for the plot e.g. "Plecoptera (HCO)". If left blank the file name is used.}
  \item{export}{Primer sequences can be exported in a text file, for further plotting. Specify a document name, if the file exists the primer data will be appended.}
  \item{rmN}{"N" nucleotides are usually not considered for coverage and base composition calculation. Set to FALS to include N Nucleotides into the calculations (See details) }
  \item{covtab}{export base poportions into file}
}
\details{
Your sequence data has to be clustered, mapped and primer regions extracted as a fasta file (With USEARCH and Geneious, see Elbrecht & Leese 2015 for details). The command PrimerPlot will only import your prepared fasta file and visualise it / calculate the degenerated primer sequence.




}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}