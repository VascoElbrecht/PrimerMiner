\name{primerversions}
\alias{primerversions}
\title{
Converting degenerated primers into their unique versions}
\description{
Degenerated primers contain variable posstions, wich makes them better matching towards different template DNAs that are commonly present in metabarcoding bulk samples. Here we take the IUPAC Codes that define the degenerated bases (e.g. \code{Y} stands for 	\code{C} or \code{T}) and convert them back to the standard bases \code{C}, \code{G}, \code{A} & \code{T}. The function takes in one string of characters (the primer sequence), and retunes all possible primer versions.

The primer sequence can also be returned as reverse complement with \code{fw=F} and iosine \code{I} bases are treated as \code{N} by default (can be turned off with \code{convertInosine=T}).

}
\usage{
primerversions(primer, fw=T, convertInosine=T)
}
\arguments{
  \item{primer}{Primer sequence in quotemarks and upper case letters.}
  \item{fw}{Primer is oriented in forward direction 5->3, should usually be the case. But also the reverse complement can be retuned by setting \code{fw=F}}
  \item{id}{
Inosine \code{I} are usually treated as \code{N} and converted into \code{C}, \code{G}, \code{A} & \code{T}. To keep the iosine, set \code{convertInosine=F}}
} 