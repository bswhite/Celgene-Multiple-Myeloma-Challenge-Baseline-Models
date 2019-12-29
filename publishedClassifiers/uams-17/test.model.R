source('uams-17.R')
source('../common/classifier-base.R')
library(GEOquery)
library(survival)

dataset <- suppressWarnings(getGEO(GEO='GSE2658', destdir=".")[[1]])
exprs(dataset) <- log2(exprs(dataset))
scores <- uams.17.eset(dataset)$res
os <- as.numeric(gsub(pattern = '.*SURTIM=(\\d+\\.?\\d*) .*',replacement = '\\1',as.character(pData(dataset)[['characteristics_ch1.2']])))
censoring <- as.numeric(gsub(pattern = '.*SURIND=(\\d) .*',replacement = '\\1',as.character(pData(dataset)[['characteristics_ch1']])))
cox <- coxph(Surv(os,event=censoring)~scores$high.risk)
print(summary(cox))
if(summary(cox)$logtest[['pvalue']] < 0.0001) {
  message('Logrank p-value matches paper')
}

data <- NULL
use.synapse <- TRUE
# Get the harmonized data
data.file <- "GSE9782APEXprobelevel_mas5.csv"
if(!file.exists(data.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7323243", downloadFile = TRUE, downloadLocation = ".")
  }
}
suppressPackageStartupMessages(library(data.table))
# data <- read.table(data.file, sep=",", header=TRUE, row.names=1)
data <- fread(data.file)
data <- as.data.frame(data)

# Drop duplicated probes
data <- data[!duplicated(data$V1),]
probe.ids <- data$V1
cols <- colnames(data)

apex.expr.data <- data[,cols[!(colnames(data) %in% c("V1"))]]
rownames(apex.expr.data) <- probe.ids


anno.file <- NULL
if(use.synapse) {
  obj <- synGet(id="syn7265181", downloadFile = TRUE, downloadLocation = ".")
  anno.file <- getFileLocation(obj)
} else {
  anno.file <- "publicClinicalHarmonized.csv"
}

anno.df <- read.table(anno.file, sep=",", header=TRUE)
rownames(anno.df) <- anno.df$SampleID
inters <- intersect(rownames(anno.df), colnames(apex.expr.data))

anno.df <- anno.df[inters,]
apex.expr.data <- apex.expr.data[,inters]

source("../../performance_metrics/metrics.R")
do.auc <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.auc(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$raw.score), digits=3)
}
do.f1 <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.f1(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$high.risk), digits=3)
}

cat("Results standardizing by test data\n")
res <- uams.17.eset(ExpressionSet(as.matrix(apex.expr.data)), standardize.using.training.data = FALSE)
classifier.res <- merge(res[["res"]], anno.df, by.x="ID", by.y="SampleID")

do.auc(classifier.res, dependent.var="PFSbasedHighRisk")
do.f1(classifier.res, dependent.var="PFSbasedHighRisk")

cat("Results standardizing by training data\n")
res <- uams.17.eset(ExpressionSet(as.matrix(apex.expr.data)), standardize.using.training.data = TRUE)
classifier.res <- merge(res[["res"]], anno.df, by.x="ID", by.y="SampleID")

do.auc(classifier.res, dependent.var="PFSbasedHighRisk")
do.f1(classifier.res, dependent.var="PFSbasedHighRisk")

trainingMeans <- c(12.7122271946531,5.70421405991235,11.2377913501713,10.1350512537729,10.2237252641476,10.1132241529945,6.83684363220325,8.4676803122732,7.89931306956123,11.3257404788118,7.66805660126567,9.07005602886575,10.9986888716969,10.5768811754997,8.62432460604888,7.81197965645517,7.04306601208372)
names(trainingMeans) <- c('200638_s_at','1557277_a_at','200850_s_at','201897_s_at','202729_s_at','203432_at','204016_at','205235_s_at','206364_at','206513_at','211576_s_at','213607_x_at','213628_at','218924_s_at','219918_s_at','220789_s_at','242488_at')
trainingSDs <- c(0.566265433826941,1.439493522169,0.497673195366107,0.934992769643516,0.935758059559669,0.684442596437629,1.22889386878789,0.904132681968134,1.01714456624349,0.91353321534587,1.53473985465217,0.942466241663803,0.485836220886129,0.760059867964268,1.75616530808258,1.36275365512788,1.42532285132829)
names(trainingSDs) <- c('200638_s_at','1557277_a_at','200850_s_at','201897_s_at','202729_s_at','203432_at','204016_at','205235_s_at','206364_at','206513_at','211576_s_at','213607_x_at','213628_at','218924_s_at','219918_s_at','220789_s_at','242488_at')

orig.trainingMeans <- trainingMeans
orig.trainingSDs <- trainingSDs

inters <- intersect(rownames(apex.expr.data), names(trainingMeans))
apex.d <- apex.expr.data[inters,]
trainingMeans <- trainingMeans[inters]
trainingSDs <- trainingSDs[inters]

testMeans <- rowMeans(apex.d)
testSDs <- apply(apex.d, 1, sd)

lm.obj <- lm(trainingMeans ~ testMeans)
r2 <- summary(lm.obj)$adj.r.squared
png("apex-vs-gse2638-means.png")
plot(testMeans, trainingMeans, main="APEX vs GSE2658 (training) means", xlab="APEX", ylab="GSE2658")
abline(lm.obj)
text(x=min(testMeans)*1.2, y=max(trainingMeans)*0.9, labels=paste0("R2 = ", signif(r2,3)))
d <- dev.off()

lm.obj <- lm(trainingSDs ~ testSDs)
r2 <- summary(lm.obj)$adj.r.squared
png("apex-vs-gse2638-sds.png")
plot(testSDs, trainingSDs, main="APEX vs GSE2658 (training) std devs", xlab="APEX", ylab="GSE2658")
abline(lm.obj)
text(x=min(testSDs)*1.2, y=max(trainingSDs)*0.9, labels=paste0("R2 = ", signif(r2,3)))
d <- dev.off()
