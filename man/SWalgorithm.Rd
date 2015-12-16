\name{SWalgorithm}
\alias{SWalgorithm}
\title{
Smith-Watermam Algorithm
}
\description{
Implement the Smith-Waterman Algorithm to produce local pairwise sequence alignment
}
\usage{
SWalgorithm(seq1, seq2, subMatrix, gapOpening = 8, gapExtension = 8)
}
\arguments{
  \item{seq1}{input the first sequence}
  \item{seq2}{input the second sequence}
  \item{subMatrix}{substitution matrix for protein sequence alignment or score matrix for DNA/RNA sequence alignment}
  \item{gapOpening}{the penalty for opening a gap in the alignment}
  \item{gapExtension}{the penalty for extending a gap by one residue}
}

\value{
  \item{path}{the optimal path in pairwaise sequence alignment}
  \item{fMatrix}{the score table used in pairwise sequence alignment to determine optimal path}
  \item{bestScore}{the best score(s) in local pairwise sequence alignemtn to start the track back}
}
\references{
Durbin, Richard, Sean R. Eddy, Anders Krogh, and Graeme Mitchison. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press, 1998.
}
\author{
Yu Gu
}