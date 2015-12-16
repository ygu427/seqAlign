\name{seqGen}
\alias{seqGen}
\title{
Generate Sequences
}
\description{
This function is used to generate multiple sequences. Type of sequence could be DNA sequence, RNA sequence and Protein sequence.
}
\usage{
seqGen(length,sampleNum,type="protein",dir="",out.file="")
}
\arguments{
  \item{length}{lenth of sequences, can be scalar or vector}
  \item{sampleNum}{number of sequences need to be generated}
  \item{type}{specify the type of sequences}
  \item{dir}{setting up the working directory for output}
  \item{out.file}{file name for output}
}
\details{
It requires package "seqinr" to write a fasta file. If user set up the working  directory and name the output file, then a fasta file contains the generated sequences would be saved. If no user-determined output file name, then there will not be a fasta file saved.
}
\value{
  \item{seqL}{the generated sequences, saved as a list object}
}

\author{
Yu Gu
Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}