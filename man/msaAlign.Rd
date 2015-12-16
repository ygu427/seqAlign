\name{msaAlign}
\alias{msaAlign}
\title{
Multiple Sequence Alignment
}
\description{
Use package "msa"
}
\usage{
msaAlign(seq, method = c("ClustalW", "ClustalOmega", "Muscle"), dir = "", out.file = "")
}
\arguments{
\item{seq}{input sequences}
\item{method}{algorithm chosen for alignment}
\item{dir}{setting up the working directory for output}
\item{out.file}{file name for output}
}
\details{
Use the package "msa"
}
\value{
\item{align}{a MsaAAMultipleAlignment object for MSA}
}

\author{
Yu Gu

Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}
\references{
Bodenhofer, Ulrich, Enrico Bonatesta, Christoph Horejs-Kainrath, and Sepp Hochreiter. "msa: an R package for multiple sequence alignment." Bioinformatics (2015): btv494.
}