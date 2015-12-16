\name{pairAlign}
\alias{pairAlign}
\title{
Pairwise Sequence Alignment
}
\description{
A common user-interface for pairwise sequence alignment. A user determined parameter input window would pop-up. User can determine using either of Needleman-Wunsch or Smith-Waterman algorithm or using both of them.
}
\usage{
pairAlign(method=c("SW","NW","both"),dir="",out.file="")
}
\arguments{
  \item{method}{algorithm chosen for alignment}
  \item{dir}{setting up the working directory for output}
  \item{out.file}{file name for output}
}

\details{
This function requires package "seqinr" and "biomaRt" and "Biostrings"
}
\value{
\item{align}{a list object contains the alignment results}
}

\author{
Yu Gu

Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}