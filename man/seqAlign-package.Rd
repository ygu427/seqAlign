\name{seqAlign-package}
\alias{seqAlign-package}
\alias{seqAlign}
\docType{package}
\title{
Sequence Alignment
}
\description{
Functions in this package implement multiple sequences alignment, such as ClustalW, ClustalOmega and Muschle, and also pairwise sequence alignbment, such as Needleman-Wuntch algorithm and Smith-waterman algorithm.
}
\details{
  \tabular{ll}{
    Package: \tab seqAlign\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2015-12-01\cr
    License: \tab GPL-3\cr
  }
To use this package, you need to update the R to at least version 3.1 and install the package "msa","Biostrings","seqinr","biomaRt","tcltk". The most important function is pairAlign and msaAlign, which are the common user interface to pairwise alignment and multiple alignment respectively.

Keywords and phrases: pairwise sequence alignment, multiple sequence alignment, Needleman-Wunsch algorithm, Smith-Waterman algorithm, Clustal W program, Clustal Omega program, MUSCLE program
}
\author{
Yu Gu

Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}
\references{
Bodenhofer, Ulrich, Enrico Bonatesta, Christoph Horejs-Kainrath, and Sepp Hochreiter. "msa: an R package for multiple sequence alignment." Bioinformatics (2015): btv494.

Durbin, Richard, Sean R. Eddy, Anders Krogh, and Graeme Mitchison. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press, 1998.
}

\keyword{ package }
