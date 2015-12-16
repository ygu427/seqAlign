msaAlign <- function(seq,method=c("ClustalW","ClustalOmega","Muscle"),
                     dir="",out.file=""){
  #require(msa)
  
  # Setting working directory and output file
  if (dir==""){
    dir <- getwd()
  } else {
    setwd(dir)
  }
  
  if (out.file==""){
    outfile <- 0
  } else outfile <- 1
  
  if (method=="ClustalW"){
    
    #######
    ## ClustalW algorithm
    #######
    align <- msa(seq,'ClustalW')
    if (outfile) {
      msaPrettyPrint(align,output='pdf',shadingMode='similar',
                     showNames='none',showLogo='none',showLegend=FALSE,
                     consensusColors='ColdHot',file=paste0(dir,'/',out.file),
                     furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                   "\\showruler{1}{top}"),
                     askForOverwrite=FALSE)
    }
  } else if (method=="ClustalOmega"){
    #######
    ## ClustalOmega algorithm
    #######
    
    align <- msa(seq,'ClustalOmega',
                               substitutionMatrix= "BLOSUM50")
    if (outfile) {
      msaPrettyPrint(align,output='pdf',shadingMode='similar',
                     showNames='none',showLogo='none',showLegend=FALSE,
                     consensusColors='ColdHot',file=paste0(dir,'/',out.file),
                     furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                   "\\showruler{1}{top}"),
                     askForOverwrite=FALSE)
    }
  } else if (method=="Muscle"){
    #######
    ## Muscle
    #######
    
    align <- msa(seq,'Muscle')
    if (outfile) {
      msaPrettyPrint(align,output='pdf',shadingMode='similar',
                     showNames='none',showLogo='none',showLegend=FALSE,
                     consensusColors='ColdHot',file=paste0(dir,'/',out.file),
                     furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                   "\\showruler{1}{top}"),
                     askForOverwrite=FALSE)
    }
  } else stop("Alignement algorithm is unrecognizable!")
  return(align)
}
