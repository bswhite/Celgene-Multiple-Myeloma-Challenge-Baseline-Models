gr_conditional_rev3 <- function(A,surv) # using A to for consistency with Prof Guan notation; using Prof Guan code as the basis
{
  t_A <- surv$t[A]
  e_A <- surv$event[A]
  r_A <- surv$r[A]
  
  if(is.na(t_A)|is.na(e_A)){return(NA)}
  surv <- na.omit(surv)
  
  cond1i <- surv$t < t_A & e_A == 1
  cond2i <- surv$t < t_A & e_A == 0
  cond3i <- surv$t > t_A & e_A == 1
  cond4i <- surv$t > t_A & e_A == 0 & surv$r==1
  cond5i <- surv$t > t_A & e_A == 0 & surv$r==0
  cond6i <- surv$t == t_A & surv$event == e_A
  cond7i <- surv$t == t_A & surv$event == 0 & e_A == 1
  
  cond         <- rep(NA,nrow(surv))
  cond[cond1i] <- 1 - (surv$r[cond1i] - r_A)/surv$r[cond1i]
  cond[cond2i] <- 0.5*(1 - (surv$r[cond2i] - r_A)/surv$r[cond2i])
  cond[cond3i] <- 1 
  cond[cond4i] <- (1 - (r_A-surv$r[cond4i])/r_A)
  cond[cond5i] <- (1 - (r_A-surv$r[cond5i])/r_A) +  0.5*(1 - (r_A-surv$r[cond5i])/r_A)
  cond[cond6i] <- 0.5
  cond[cond7i] <- 1
  
  return(sum(cond,na.rm=T))
}

guan_rank <- function(survX) # Guan_rank transform  takes in a Surv object and transforms it to the Prof Yuanfang Guan's full rank as used in ALS patient stratification challenge
{
  suppressPackageStartupMessages(library("survival"))
  kmFit  <- survfit(survX~1)
  survDF <- data.frame(t=survX[,1], event = survX[,2], r = kmFit$surv[match(survX[,1],kmFit$time)])
  gRank  <- sapply(1:nrow(survDF), gr_conditional_rev3, surv=survDF)
  
  return(gRank)
}


