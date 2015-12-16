### SWalgorithm.R
###
### Implement the dynamic programming algorithm
### Smith-Waterman algorithm for local alignment
### Pairwise protein sequence local alignment
### User input: two GeneBank identifier for proteins
###          or input sequences directly (check box)
###          or read FASTA file (if no input for identifieer)
###          gap-opening and gap-extension penalties
###          score matrix 
### 
### Result: local optimal paths (list for multiple paths)
###         the best score
###         the F Score Matrix
### Result would be written into a text file
###
### Required packages:
###         "biomaRt" -- for sequence query
###         "tcltk" -- for user input interface
###         "seqinr" -- for read FASTA file 
###
### Written by Yu Gu ; Supervised by Hongyu Miao
### 03-17-2015


###############################################
##                                           ##
##       Local Alignment Algorithm           ##
##                                           ##
###############################################

SWalgorithm <- function(seq1,seq2,subMatrix,gapOpening = 8, gapExtension = 8){
  
  ###
  #
  # Conver the sequences to upper case. 
  # Doing this to make sure all chars in the string is upper case
  #
  ###
  toupper(seq1)
  toupper(seq2)
  
  ###
  #
  # Create vectors and matrices to store temp and final results
  # alignX & alignY store sequences as character vector
  # fMatrix: stores the scores
  # track: stores the postion from which getting to the best score
  #        index on fMatrix
  # gapEx & gapEy: store number of gap extension 
  #
  ###
  lenX <- nchar(seq1)
  lenY <- nchar(seq2)
  
  alignX <- vector()
  alignY <- vector()
  
  for (i in 1:lenX) {
    alignX[i] <- substr(seq1,i,i)
  }
  
  for (j in 1:lenY) {
    alignY[j] <- substr(seq2,j,j)
  }
  
  fMatrix <- matrix(0,lenY+1,lenX+1)
  
  track <- list()
  gapEx <- matrix(0,lenY,lenX)
  gapEy <- matrix(0,lenY,lenX)
  
  ###
  #
  # Initialization:
  # initial coordinate (1,1), but for fMatrix, it's F(2,2)
  # bestScore = 0 indicates starting here
  #
  ###
  
  bestScore <- 0
  
  ###
  #
  # Filling out the F Matrix
  # Note the current cell is F(j+1,i+1)
  # scoreI: F(j,i+1) -> F(j+1,i+1)
  # scoreJ: F(j+1,i) -> F(j+1,i+1)
  # F(1,) = F(0,1) = 0
  # if F(j,i)<0, set F(j,i)=0
  #
  ###
  for (j in 1:lenY) {
    for (i in 1:lenX){
      
      # get the symbols from both sequences 
      X <- alignX[i]
      Y <- alignY[j]
      
      # score from diagonal
      s <- subMatrix[X,Y]
      score <- fMatrix[j,i]+s
      if (score<0) score<-0
      bestScore <- score
      pair <- as.character(c(j,i))
      
      # score from upper cell, which causes gap in Xi
      # gapX stores the number of gap extension on X
      # compared with current best score from disgonal:
      #   if better, replace the current one with the new score
      #   if same, keep both of the current one and the new one
      #   if worse, discard the new score
      
      if (j==1) {
        gapX <- 0
      } else {
        gapX <- gapEx[j-1,i]
      }
      if (gapX==0) {
        scoreI <- fMatrix[j,i+1] - gapOpening
      } else {
        scoreI <- fMatrix[j,i+1] - gapX*gapExtension
      }
      if (scoreI<0) scoreI <- 0
      if (scoreI > bestScore) {
        bestScore <- scoreI
        pair <- c(j,i+1)
        gapEx[j,i] <- gapX+1
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-"X"
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } else if (scoreI==bestScore) {
        pair <-c(pair,c(j,i+1))
        gapEx[j,i] <- gapX+1
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-c(NA,"X")
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } else {
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-NA
        gapEx[j,i]<-0
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } # end of checking upper cell
      
      # score from left cell, which causes gap in Yj
      # gapY store the number of gap extension on Y
      # compared with the current best score (diagonal, upper, or both)
      # similar process as the one for upper cell
      if (i==1) {
        gapY <- 0
      } else {
        gapY <- gapEy[j,i-1]
      }
      if (gapY==0) {
        scoreJ <- fMatrix[j+1,i] - gapOpening
      } else {
        scoreJ <- fMatrix[j+1,i] - gapY*gapExtension
      }
      if (scoreJ<0) scoreJ <- 0
      if (scoreJ > bestScore) {
        bestScore <- scoreJ
        pair <- c(j+1,i)
        gapEy[j,i] <- gapY+1
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-"Y"
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } else if (scoreJ==bestScore) {
        pair <-c(pair,c(j+1,i))
        gapEy[j,i] <- gapY+1
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-c(gap,"Y")
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } else {
        map <- matrix(pair,ncol=2,byrow=TRUE)
        colnames(map)<-c("J","I")
        gap<-gap
        gapEy[j,i]<-0
        map<-cbind(map,gap)
        track[[(j-1)*lenX+i]]<-map
      } # end of checking left cell
      
      fMatrix[j+1,i+1]<-max(bestScore,0)
      #track
    } # end of i
  } # end of j
  
  ###
  #
  # Track back
  # done <logic value>: 0 indicates track back not done yet; 
  #                     1 indicates track back done 
  # get the maximum score(s) in fMatrix
  # for each maximum score, track back the path indexes
  # retrieve gap info. from track list
  # fill out the path matrix
  #
  ###
  
  fScore <- max(fMatrix)
  start.matrix <- which(fMatrix==fScore,arr.ind=TRUE)
  nStart <- nrow(start.matrix)
  path.list <- list()
  
  for (init in 1:nStart) {
    
    # Initial Status
    done <- 0
    start.index <- start.matrix[init,]
    x.index <- start.index[2]
    y.index <- start.index[1]
    path.index <- matrix(c(y.index,x.index),ncol=2)
    count <- as.integer((y.index-1-1)*lenX + x.index-1)
    
    # find the local optimal path
    while(!done) {
      trBack <- track[[count]]
      x.index <-as.integer(trBack[2])
      y.index <-as.integer(trBack[1])
      sr <- fMatrix[y.index,x.index]
      if (sr==0){
        done <- 1
        break
      }
      path.index <- rbind(path.index,c(y.index,x.index))
      count <- (y.index-1-1)*lenX + x.index-1 
    } # end of while
    l<-nrow(path.index)
    for (k in 1:floor(l/2)) {
      temp <- path.index[k,]
      path.index[k,] <- path.index[(l-k+1),]
      path.index[(l-k+1),] <- temp
    }
    
    path <- matrix(NA,ncol=2)
    colnames(path)<-c("X.align","Y.align")
    
    for (m in 1:l) {
      symX.index <-path.index[m,2]
      symY.index <-path.index[m,1]
      count <- (symY.index-1-1)*lenX + symX.index-1 
      gapD <- track[[count]][3]
      if (is.na(gapD)) {
        X <- alignX[symX.index-1]
        Y <- alignY[symY.index-1]
        path <- rbind(path,c(X,Y))
      } else if (gapD == 'X') {
        X <- '_'
        Y <- alignY[symY.index-1]
        path <- rbind(path,c(X,Y))
      } else {
        X <- alignX[symX.index-1]
        Y <- '-'
        path <- rbind(path,c(X,Y))
      }
    }
    path <- path[-1,]
    path.list[[init]] <- t(path)
  } # end of for
  
  rownames(fMatrix) <- c(" ", alignY)
  colnames(fMatrix) <- c(" ", alignX)
  return (list(path.list = path.list, bestScore = fScore, fMatrix = fMatrix))
}