library(RUnit)
library(seqAlign)

path <- "C:/Users/ygu/Documents/GitHub/seqAlign/R"
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

dir<-sourceDir(path)
test.suite <- defineTestSuite("example",
                              dirs = file.path("C:/Users/ygu/Documents/GitHub/seqAlign/tests"),
                              testFuncRegexp = "^test.+",
                              testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)