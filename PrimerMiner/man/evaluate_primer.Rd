\name{evaluate_primer}
\alias{evaluate_primer}
\title{
in silico primer evaluation against sequences
}
\description{
Compare primer sequence against each sequence in a sequence alignment (supply fasta file). Primer sequence, binding site and table giving possition dependent penalty values have to be supplied. Adjacent nucleotides and mismatch types (e.g. A-A or A-C etc) is also considdered.

For each sequence penalty scores are calculated and saved in a csv file. Penalty scores can be further analysed statisticly or all sequences above a certain treshold be considdered not working. See \link{primer_threshold} for details.
}
\usage{
evaluate_primer(alignment_imp, primer_sequ, start, stop, forward = T, save = NULL, gap_NA = T, N_NA = T, mm_position = NULL, mm_type = NULL, adjacent = 2, sequ_names=T)
}
\arguments{
  \item{alignment_imp}{
Sequence alignment as a fasta file, which the primer should be evaluated against.
}
  \item{primer_sequ}{
Primer sequence (5 -> 3\'), can contain wobble bases and inosine (I).
}
  \item{start}{
Starting possition of the primer binding region in the alignment.
}
  \item{stop}{
End possition of the primer binding region in the alignment.
}
  \item{forward}{
Indicate if the primer binds forward (default = T), or in reverse complement (set to F).
}
  \item{save}{
Save the result table under this file name (csv).
}
  \item{gap_NA}{
If TRUE gaps (-) are replaced with NA for missing data, if FALSE gaps will be counted as mismatches and penalty scores will be applied. 
}
  \item{N_NA}{
If TRUE ambigious bases (N) are replaced with NA for missing data, if FALSE gaps will be counted as mismatches and penalty scores will be applied.
}
  \item{mm_position}{
Suppliy a table here indicating mismatch penalty scores for each possition (Primer tip 3’ to 5’ tail). See https://github.com/VascoElbrecht/PrimerMiner/tree/master/primer_scoring for examples.
}
  \item{mm_type}{
Suppliy a table here indicating penalty adjustments depending on the mismatch type (scores get multiplied with the respective value). See https://github.com/VascoElbrecht/PrimerMiner/tree/master/primer_scoring for examples.
}
  \item{adjacent}{
Factor with wich penalty scores of adjacent nucleotides should be multiplied. To turn this option off set to 1.
}
  \item{sequ_names}{
If true, includes template sequence names as provided in the sequence alignment (fasta file).
}
}
\details{
}
\value{
}
\references{
Elbrecht & Leese (2016). PrimerMiner: an r package for development and in silico validation of DNA metabarcoding primers. Methods in Ecology and Evolution \href{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12687/full}{DOI: 10.1111/2041-210X.12687}
}
\author{
}
\note{
}

\seealso{
\code{\link{primer_threshold}}
}
\examples{
}
