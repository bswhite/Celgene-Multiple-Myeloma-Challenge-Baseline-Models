
iss.stage <- function(X, y, threshold=2) {
  iss.to.numeric <- function(X) {
    convert <- c(1,2,3);
    names(convert) <- c('I','II','III')
    if(!is.na(X)) {
      convert[X]
    } else {
      NA
    }
  }
  iss <- unlist(lapply(X$ISSstage,iss.to.numeric))
  res <- data.frame(ID=sampleNames(test.eset),raw.score=iss,high.risk=(iss > threshold))
  obj <- list(threshold=threshold)
  return(list(obj=obj,res=res))
}
