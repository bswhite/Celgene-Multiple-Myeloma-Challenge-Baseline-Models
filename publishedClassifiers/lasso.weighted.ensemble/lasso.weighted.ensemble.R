lasso.weighted.ensemble <- function(X, y, threshold=NULL) {

  emc92res <- emc.92.entrez(X = X, y = NULL, threshold = NULL)$res
  uams5res <- uams.5.entrez(X = X, y = NULL, threshold = NULL)$res
  uams17res <- uams.17.entrez(X = X, y = NULL, threshold = NULL)$res
  uams70res <- uams.70.entrez(X = X, y = NULL, threshold = NULL)$res
  
  raw.score <- 0.0391*u70res$raw.score+0.1959*u17res$raw.score+-0.1168*u5res$raw.score+0.5091*emc92res$raw.score
  if(is.null(threshold)) {
    threshold = median(raw.score)
  }
  df <- data.frame(ID=emc92res$ID,raw.score=raw.score,high.risk=round(raw.score,digits=2) > threshold)
  obj <- list(threshold=threshold)
  return(list(obj=obj,res=df))
}

lasso.weighted.ensemble.iss <- function(X, y, threshold=NULL) {
  res <- lasso.weighted.ensemble(X[,!(colnames(X) %in% "ISSstage")], y, threshold)$res

  if(is.null(threshold)) {
    threshold = median(res$raw.score)
  }
  df <- data.frame(ID=res$ID,raw.score=res$raw.score,high.risk=round(res$raw.score,digits=2) > threshold)
  df$high.risk <- df$high.risk | X$ISSstage %in% c('III', 'II')
  obj <- list(threshold=threshold)
  return(list(obj=obj,res=df))
}