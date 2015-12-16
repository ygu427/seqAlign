### pairAlign.R
### pairwise sequence alignment for protein sequences
### input: the type of methods
###           "SW"--Smith-Waterman Algorithm for local alignment
###           "NW"--Needleman-Wunsch Algorithm for global alignment
###           "both"--will run both of above algorithms
###        dir -- for working directory
###        out.file -- output file, if NA, then no output will be saved
### output: alignment results
###         if out.file is input, then the output would be saved under
###         the user-defined directory
### Written by Yu Gu
### Dec. 6th, 2015

pairAlign <- function(method=c("SW","NW","both"),dir="",out.file=""){
  
  #require("seqinr")
  #require("biomaRt")
  #require("Biostrings")
  
  ### Checking input parameters
  if (method!="SW" | method!="NW" | method!="both"){
    stop("Alignement algorithm is unrecognizable!")
  }
  
  # Setting working directory and output file
  if (dir==""){
    dir <- getwd()
  } else {
    setwd(dir)
  }
  
  if (out.file==""){
    outfile <- 0
  } else outfile <- 1
  
  input <- varEntryDialog()
  
  if (input$seqCB==1) {
    seq1 <- input$Sequence1
    seq2 <- input$Sequence2
  } else if (input$FASTA==1) {
    #seq1 <- read.fasta(input$Identifier1,forceDNAtolower = FALSE)
    #seq2 <- read.fasta(input$Identifier2,forceDNAtolower = FALSE)
    infile <- read.fasta(input$Sequence1,seqtype = "AA",seqonly = TRUE)
    seq1 <- infile[[1]]
    seq2 <- infile[[2]]
  } else {
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    seq1 <- getSequence(id=input$Sequence1,type="embl",
                        seqType="peptide", mart=ensembl)
    seq1 <- seq1$peptide
    seq2 <- getSequence(id=input$Sequence2,type="embl",
                        seqType="peptide", mart=ensembl)
    seq2 <- seq2$peptide
  }
  
  gapOpening <- input$gapPenalty
  gapExtension <- input$extensionPenalty
  substitutionMatrix <- input$subMatrix
  
  data(substitutionMatrix)
  if (substitutionMatrix=="BLOSUM62") {
    subMatrix <- BLOSUM62
  } else if (substitutionMatrix=="BLOSUM50") {
    subMatrix <- BLOSUM50
  } else if (substitutionMatrix=="PAM250") {
    subMatrix <- PAM250
  } else if (substitutionMatrix=="BLOSUM45") {
    subMatrix <- BLOSUM45
  } else if (substitutionMatrix=="BLOSUM80") {
    subMatrix <- BLOSUM80
  } else if (substitutionMatrix=="PAM120") {
    subMatrix <- PAM120
  }
  
  align <- list()
  
  if (method=="SW"){
    swAlign <- SWalgorithm(seq1,seq2,subMatrix,gapOpening,gapExtension)
    return(align=swAlign)
  } else if(method=="NW"){
    nwAlign <- NWalgorithm(seq1,seq2,subMatrix,gapOpening,gapExtension)
    return(align=nwAlign)
  } else {
    swAlign <- SWalgorithm(seq1,seq2,subMatrix,gapOpening,gapExtension)
    nwAlign <- NWalgorithm(seq1,seq2,subMatrix,gapOpening,gapExtension)
    return(align=list(swAlign=swAlign,nwAlign=nwAlign))
  }
  
  if (outfile) {
    sink(paste0(dir,'/',out.file))
    output <- pairAlign()
    print(output)
    sink()
  }
} 