### seqGen.R
###
### Applying this function to generate sample sequences
### Could be DNA or RNA or Protein sequences
### Input: 
###   lenth -- lenth of sequences, can be scalar or vector
###   sampleNum -- number of sequences need to be generated
###   type -- specify the type of sequences
###   dir -- set up the working directory
###   out.file -- output file's name
### Output:
###   generated sequences
###   true alignment
###   fasta file saved under the user-defined directory or
###   current working directory if no user-defined one
###
### Written by Yu Gu
### Dec. 1st 2015

seqGen <- function(length,sampleNum,type="protein",dir="",out.file=""){
  #require("seqinr")
  
  if (dir==""){
    dir <- getwd()
  } else {
    setwd(dir)
  }
  
  if (out.file==""){
    outfile <- 0
  } else outfile <- 1
  
  ## Letter repository for sequences
  DNA <- c("A","T","C","G")
  RNA <- c("A","U","C","G")
  amino <- c("A","B","C","D","E","F","G","H","I","K","L","M",
             "N","P","Q","R","S","T","V","W","X","Y","Z")
  tank <- vector()
  if (type=="dna"){
    tank <- DNA
  } else if (type=="rna"){
    tank <- RNA
  } else {
    tank <- amino
  }
  
  original.seq <- sample(tank,length,T)
  seqS <- matrix(original.seq,sampleNum,length,T)
  seqL <- list()
  seqL[[1]] <- original.seq
  diffM <- matrix(0,sampleNum,length)
  trueAlign <- vector()
  for (i in 2:sampleNum) {
    temp <- seqS[i,]
    startP <- ceiling(runif(1,0,5))
    ind <- startP
    repeat {
      interval <- sample(1:3,1,T,prob=c(0.85,0.1,0.05))
      if ((ind+interval-1) >= length) break()
      replace <- sample(tank,interval,T)
      position <- seq(from=ind,by=1,length.out=interval)
      temp[position] <- replace
      diffM[i,position] <- 1
      ind <- ind + interval + ceiling(runif(1,20,30))
      if (ind >= length) break()
    }
    seqS[i,] <- temp
    seqL[[i]] <- temp
  }
  index <- which(colSums(diffM)==0)
  trueAlign[index] <- original.seq[index]
  trueAlign[-index] <- "-"
  if (outfile) {
    write.fasta(seqL,names=paste0("Seq",seq(1:sampleNum)),
                file.out=out.file)
  }
  return(seqL)
}