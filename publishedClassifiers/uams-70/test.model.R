source('uams-70.R')
source('../common/classifier-base.R')
library(GEOquery)
library(survival)

dataset <- suppressWarnings(getGEO(GEO='GSE57317')[[1]])
exprs(dataset) <- log2(exprs(dataset))
u70.scores <- uams.70.eset(dataset)$res
gep70.score <- as.numeric(gsub(pattern = 'gep70 score: ',replacement = '',as.character(dataset[['characteristics_ch1.5']])))
if(sum(gep70.score-round(u70.scores$raw.score,digits=4)) < 10^(-4)) {
  message('GEP70 scores match paper')
}

