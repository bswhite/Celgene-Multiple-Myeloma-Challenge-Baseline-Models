naive.ensemble <- function(X, y, threshold=NULL) {

  emc92res <- emc.92.entrez(X = X, y = NULL, threshold = NULL)$res
  u5res <- uams.5.entrez(X = X, y = NULL, threshold = NULL)$res
  u17res <- uams.17.entrez(X = X, y = NULL, threshold = NULL)$res
  u70res <- uams.70.entrez(X = X, y = NULL, threshold = NULL)$res
  
  raw.score <- u70res$raw.score*u17res$raw.score*u5res$raw.score*emc92res$raw.score
  if(is.null(threshold)) {
    threshold = median(raw.score)
  }
  df <- data.frame(ID=emc92res$ID,raw.score=raw.score,high.risk=round(raw.score,digits=2) > threshold)
  obj <- list(threshold=threshold)
  return(list(obj=obj,res=df))
}
