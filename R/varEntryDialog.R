###############################################
##                                           ##
##             User Input GUI                ##
##                                           ##
###############################################

###
#
# Return: GeneBank Identifiers or Sequences or FASTA file directory
#         gap openting & extension penalty
#         Substitution Matrix
#
###
varEntryDialog <- function() {
  #require("tcltk")
  
  ###
  #
  # variables and title
  # vars <- c("Sequence1","Sequence2","gap opening","gap extension",
  #          "Substitution Matrix","check sequence","check FASTA file")
  #
  ###
  
  title <- "input parameters"
  
  ###
  #
  # Create a variable to keep track of the state of the dialog window:
  # done = 0; If the window is active
  # done = 1; If the window has been closed using the OK button
  # done = 2; If the window has been closed using the Cancel button or destroyed
  #
  ###
  done <- tclVar(0)
  
  tt <- tktoplevel()
  tkwm.title(tt, title)  
  entries <- list()
  tclvars <- list()
  
  tkbind(tt,"<Destroy>",function() tclvalue(done)<-2)
  
  ###
  #
  # Initial values
  #
  ###
  tclvars[[1]] <- tclVar("")
  entries[[1]] <- tkentry(tt, textvariable=tclvars[[1]])
  
  tclvars[[2]] <- tclVar("")
  entries[[2]] <- tkentry(tt, textvariable=tclvars[[2]])
  
  tclvars[[3]] <- tclVar("8")
  entries[[3]] <- tkentry(tt, textvariable=tclvars[[3]])
  
  tclvars[[4]] <- tclVar("8")
  entries[[4]] <- tkentry(tt, textvariable=tclvars[[4]])
  
  ###
  #
  # Check boxes: whether identifier inputs are Sequences
  #
  ###
  
  seq.CB <- tkcheckbutton(tt)
  seq.cbValue <- tclVar("0")
  tkconfigure(seq.CB,variable=seq.cbValue)
  
  ###
  #
  # Choose Substitution Matrix
  # Default one is BLOSUM50.  Indexing starts at zero.
  #
  ###
  
  scr <- tkscrollbar(tt, repeatinterval=5,
                     command=function(...) tkyview(tl,...))
  tl<-tklistbox(tt,height=4,selectmode="single",
                yscrollcommand=function(...) tkset(scr,...),
                background="white")
  subMatrix <- c("BLOSUM45","BLOSUM50","BLOSUM62","BLOSUM80","PAM120","PAM250")
  for (i in (1:6))
  {
    tkinsert(tl,"end",subMatrix[i])
  }
  tkselection.set(tl,1)
  
  ###
  #
  # Process : submit & cancel & reset
  #
  ##
  
  doneVal <- as.integer(tclvalue(done))
  results <- list()
  
  # RESET
  reset <- function() {
    tclvalue(tclvars[[1]]) <<- ""
    tclvalue(tclvars[[2]]) <<- ""
    tclvalue(tclvars[[3]]) <<- "8"
    tclvalue(tclvars[[4]]) <<- "8"
    seq.cbValue <- tclVar("0")
    tkconfigure(seq.CB,variable=seq.cbValue)
  }
  reset.but <- tkbutton(tt, text="Reset", command=reset)
  
  # CANCEL
  cancel <- function() {
    tclvalue(done) <- 2
  }
  cancel.but <- tkbutton(tt, text='Cancel', command=cancel)
  
  # SUBMIT
  submit <- function() {
    tryCatch( {
      results[["Sequence1"]] <<- tclvalue(tclvars[[1]])
      results[["Sequence2"]] <<- tclvalue(tclvars[[2]])
      results[["gapPenalty"]] <<- as.integer(tclvalue(tclvars[[3]]))
      results[["extensionPenalty"]] <<- as.integer(tclvalue(tclvars[[4]]))
      cbVal <- as.integer(tclvalue(seq.cbValue))
      results[["seqCB"]] <<- cbVal
      matrixChoice <- subMatrix[as.integer(tkcurselection(tl))+1]
      results[["subMatrix"]] <<- matrixChoice
      results[["FASTA"]] <<- 0
      tclvalue(done) <- 1
    },
    error = function(e) { tkmessageBox(message=geterrmessage()) },
    finally = { }
    )
  }
  submit.but <- tkbutton(tt, text="Submit", command=submit)
  
  
  ###
  #
  # construct the GUI
  #
  ###
  
  # parameter input
  heading1 <- tklabel(tt,text="Input Sequences")
  tkgrid(heading1,columnspan=2,pady=10)
  tkgrid.configure(heading1,sticky='we')
  seq1 <- tklabel(tt,text="Sequence 1")
  tkgrid(seq1,entries[[1]],pady=10, padx=10)
  seq2 <- tklabel(tt,text="Sequence 2")
  tkgrid(seq2,entries[[2]],pady=10, padx=10)
  is.seq <- tklabel(tt,text="Sequence?")
  tkgrid(is.seq,seq.CB)
  tkgrid.configure(seq1,seq2,is.seq,sticky='e')
  tkgrid.configure(entries[[1]],entries[[2]],seq.CB,sticky='w')
  
  heading2 <- tklabel(tt,text="Please enter the penalties")
  tkgrid(heading2,columnspan=2,pady=10)
  tkgrid.configure(heading2,sticky='we')
  gap <- tklabel(tt,text="gap")
  extension <- tklabel(tt,text="extension")
  tkgrid(gap,entries[[3]],pady=10, padx=4)
  tkgrid(extension,entries[[4]],pady=10, padx=4)
  tkgrid.configure(gap,extension,sticky='e')
  tkgrid.configure(entries[[3]],entries[[4]],sticky='w')
  
  heading3 <- tklabel(tt,text="Choose the Substitution Matrix")
  tkgrid(heading3,columnspan=2,pady=10)
  tkgrid.configure(heading3,sticky='we')
  tkgrid(tl,scr)
  tkgrid.configure(scr,rowspan=4,sticky='nsw')
  
  # Processing commands
  tkgrid(submit.but, cancel.but, reset.but, pady=10, padx=10)
  tkfocus(tt)
  
  # Do not proceed with the following code until the variable done is non-zero.
  #   (But other processes can still run, i.e. the system is not frozen.)
  tkwait.variable(done)
  
  if(tclvalue(done) != 1) {
    results <- NULL
  }
  
  tkdestroy(tt)
  if (results[["Sequence1"]]=="" | results[["Sequence2"]]=="") {
    tkmessageBox(message = "Please choose the first FASTA file")
    fileName <- tclvalue(tkgetOpenFile()) 
    if (!nchar(fileName)) {
      tkmessageBox(message = "No file was selected!")
    } else {
      tkmessageBox(message = paste("The file selected was", fileName))
    }
    results[["Sequence1"]] <- fileName
    results[["FASTA"]] <- 1
  }
  
  #if (results[["Sequence2"]]==""){
  #tkmessageBox(message = "Please choose the second FASTA file")
  #fileName <- tclvalue(tkgetOpenFile()) # Very simple, isn't it?
  #if (!nchar(fileName)) {
  #tkmessageBox(message = "No file was selected!")
  #} else {
  #tkmessageBox(message = paste("The file selected was", fileName))
  #}
  #results[["Sequence2"]] <- fileName
  #results[["FASTA"]] <- 1
  #}
  return(results)
}