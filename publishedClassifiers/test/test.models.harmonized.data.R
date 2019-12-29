source("../common/classifier-base.R")
source("../emc-92/emc-92.R")
source("../uams-5/uams-5.R")
source("../uams-17/uams-17.R")
source("../uams-70/uams-70.R")
source("../iss-stage/iss.stage.R")

use.synapse <- TRUE

seed <- 1234
# seed <- 2334
cat(paste0("seed: ", seed, "\n"))
set.seed(seed)

if(use.synapse) {
  library(synapseClient)
  synapseLogin()
}

# This follows http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# counts has genes as rows and samples as columns
get.variable.genes <- function(counts) {
  library(DESeq); library(statmod); library(pcaMethods); library(fastICA)
  require(DESeq)
  lib.size <- estimateSizeFactorsForMatrix(counts)
  ed <- t(t(counts)/lib.size)
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  pdf("variable.pdf")
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
  smoothScatter(log(means),log(cv2))
  
  require(statmod)
  minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
  minMeanForFit <- 0.001
  useForFit <- means >= minMeanForFit # & spikeins
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  
  # repeat previous plot
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
  xg <- exp(seq( min(log(means[means>0])), max(log(means[means>0])), length.out=1000 ))
  vfit <- a1/xg + a0
  # add fit line
  lines( log(xg), log(vfit), col="black", lwd=3 )
  df <- ncol(ed) - 1
  # add confidence interval
  lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
  lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")  
  
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio,decreasing=T)
  oed <- ed[varorder,]

  # repeat previous plot
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
  # add top 100 genes

  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  sigVariedGenes <- adj.pval<1e-3;
  points(log(means[sigVariedGenes]),log(cv2[sigVariedGenes]),col=2)
  # print(table(sigVariedGenes))  
  
  d <- dev.off()
  return(names(means)[sigVariedGenes])
}

data <- NULL
# Get the harmonized data
data.file <- "studyCorrectedPubicData.csv"
if(!file.exists(data.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7517756", downloadFile = TRUE, downloadLocation = ".")
  }
}

suppressPackageStartupMessages(library(data.table))
data <- fread(data.file)
data <- as.data.frame(data)
entrez.ids <- data$V1

entrez.ids <- data$V1
cols <- colnames(data)
expr.data <- data[,cols[!(colnames(data) %in% c("V1"))]]
rownames(expr.data) <- entrez.ids

# Some of the columns seem to be rederived from the same patient, with names like
# MMRF_2087_1_BM and MMRF_2087_2_BM
# or
# MMRF_1992_1_BM and MMRF_1992_1_PB

# Let's take a single column per patient, favoring "*_1_BM" over "*_1_PB" over "*_2_BM"
# Do this by sorting the columns by name and taking the first non-duplicated.
expr.data <- expr.data[,order(colnames(expr.data))]
orig.cols <- colnames(expr.data)
cols <- orig.cols
cols <- gsub(cols, pattern="_\\d+_BM", replacement="")
cols <- gsub(cols, pattern="_\\d+_PB", replacement="")

## cat(paste0("Origin matrix had ", ncol(expr.data), " cols\n"))
expr.data <- expr.data[,!duplicated(cols, fromLast=FALSE)]
orig.cols.keep <- orig.cols[!duplicated(cols, fromLast=FALSE)]
expr.samples <- data.frame(orig=orig.cols.keep, sanitized=cols)
## cat(paste0("After dropping duplicate samples, matrix has ", ncol(expr.data), " cols\n"))
        
orig.expr.cols <- colnames(expr.data)
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_BM", replacement="")
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_PB", replacement="")
expr.cols <- colnames(expr.data)

# Pull in the MMRF annotation data.  
anno.file <- NULL
if(use.synapse) {
  obj <- synGet(id="syn7488088", downloadFile = TRUE, downloadLocation = ".")
  anno.file <- getFileLocation(obj)
} else {
  anno.file <- "MMRF patient-level clinical and cytogenetic inventory Harmonized.csv"
}

anno.df <- read.table(anno.file, sep=",", header=TRUE)
anno.df$ID <- apply(anno.df[,c("Study","SampleID")], 1, function(row) paste(row[1], row[2], sep=" "))
# This must have columns event, time, and ID.  
# Use PFS annotations to define high risk by setting event = Progression; time = PFStimeMonths; ID = SampleID
anno.df <- anno.df[,c("ID", "OStimeMonths", "OSdeath", "PFStimeMonths", "Progression", "OSbasedHighRisk", "PFSbasedHighRisk", "Sex", "ISSstage", "Age")]
colnames(anno.df) <- c("ID", "OStime", "OSevent", "time", "event", "OSbasedHighRisk", "PFSbasedHighRisk", "Sex", "ISSstage", "Age")

# Pull in the public clinical data
public.anno.file <- NULL
if(use.synapse) {
  obj <- synGet(id="syn7265181", downloadFile = TRUE, downloadLocation = ".")
  public.anno.file <- getFileLocation(obj)
} else {
  public.anno.file <- "publicClinicalHarmonized.csv"
}

public.anno.df <- read.table(public.anno.file, sep=",", header=TRUE)
public.anno.df$ID <- apply(public.anno.df[,c("Study","SampleID")], 1, function(row) paste(row[1], row[2], sep=" "))
# This must have columns event, time, and ID.  
# Use PFS annotations to define high risk by setting event = Progression; time = PFStimeMonths; ID = SampleID
public.anno.df <- public.anno.df[,c("ID", "OStimeMonths", "OSdeath", "PFStimeMonths", "Progression", "OSbasedHighRisk", "PFSbasedHighRisk", "Sex", "ISSstage", "Age")]
colnames(public.anno.df) <- c("ID", "OStime", "OSevent", "time", "event", "OSbasedHighRisk", "PFSbasedHighRisk", "Sex", "ISSstage", "Age")

if(!all(colnames(public.anno.df) == colnames(anno.df))) {
  cat("Column names do not match\n")
}

anno.df <- rbind(anno.df, public.anno.df)
anno.df$orig.ID <- anno.df$ID
anno.df$ID <- unlist(lapply(anno.df$ID, function(x) gsub(x=x, pattern="GSE19784HOVON65", replacement="HOVON65")))
anno.df$ID <- unlist(lapply(anno.df$ID, function(x) gsub(x=x, pattern="-", replacement=".")))

flag <- colnames(expr.data) %in% anno.df$ID
# colnames(expr.data)[!flag]

## HERE

# Exclude any samples that do have a PFS-based high risk
## cat(paste0("Excluding ", length(which(is.na(anno.df$PFSbasedHighRisk))), " of ", nrow(anno.df), " samples without low/high risk\n"))
anno.df <- anno.df[!is.na(anno.df$PFSbasedHighRisk),]

orig.anno.samples <- nrow(anno.df)
orig.expr.samples <- ncol(expr.data)
rownames(anno.df) <- anno.df$ID
orig.row.names <- rownames(anno.df)
orig.expr.names <- colnames(expr.data)
inters <- intersect(rownames(anno.df), colnames(expr.data))
cat(paste0(length(inters), " columns are shared between the ", nrow(anno.df), " clinically annotated samples and the ", ncol(expr.data), " expression samples\n"))

colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_BM", replacement="")
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_PB", replacement="")
expr.cols <- colnames(expr.data)


# Subset the data and annotations to include only samples appearing in both.
anno.df <- anno.df[inters,]
expr.data <- expr.data[,inters]

## Filter non-expressed genes (unless they are involved in a classifer):
coefficients <- c(emc.92.get.coefficients(), uams.17.get.coefficients(), uams.70.get.coefficients(), uams.5.get.coefficients())

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
probesets <- names(coefficients)
bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), 
            filters = 'affy_hg_u133_plus_2', 
            values = probesets, 
            mart = ensembl)
names(bm) <- c("PROBE", "GENE")
coefficient.genes <- bm$GENE

anno.df$ID <- as.character(anno.df$ID)
expr.tbl <- expr.data[,anno.df$ID]
eset <- ExpressionSet(as.matrix(expr.tbl))

translocation.targets <- c("FGFR3", "MMSET", "WHSC1", "CCND3", "CCND1", "MAF", "MAFB")

highly.mutated.genes <- c("CDKN2C", "RB1", "CCND1", "CDKN2A", "NRAS", "KRAS", "BRAF", "MYC", "TRAF3", "CYLD", "DKK1", "FRZB", "DNAH5", "XBP1", "BLIMP1", "PRDM1", "IRF4", "TP53", "MRE11A", "PARP1", "DIS3", "FAM46C", "LRRK2", "KDM6A", "UTX", "MLL", "MMSET", "WHSC1", "HOXA9", "KDM6B", "IGLL5")

sym.to.id <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), 
                       filters = 'hgnc_symbol', 
                       values = c(translocation.targets, highly.mutated.genes), 
                       mart = ensembl)
names(sym.to.id) <- c("SYMBOL", "GENE")

## When looking at the expression data, check genes that are highly mutated/translocated or are used by one
## of the other classifiers, and are expressed
expr.gene.subset <- intersect(rownames(eset), sym.to.id$GENE)
expr.gene.subset <- c(expr.gene.subset, coefficient.genes)
## expr.gene.subset <- c(expr.gene.subset, get.variable.genes(expr.data))
expr.gene.subset <- intersect(expr.gene.subset, rownames(eset))
cat(paste0("Number of expression features: ", length(expr.gene.subset), "\n"))

## Restrict to highly variable genes
## Use every expression feature
## expr.gene.subset <- rownames(eset)

# Finally, make sure we don't have any NAs in any of our features.  If we do, drop the corresponding sample.
na.samples <- c()
flag <- unlist(apply(anno.df[,c("Sex", "Age", "ISSstage")], 1, function(row) any(is.na(row)) || any(row=="")))
if(any(flag)) {
  na.samples <- c(na.samples, rownames(anno.df)[flag])
}

flag <- unlist(apply(exprs(eset)[expr.gene.subset,], 2, function(col) any(is.na(col))))
if(any(flag)) {
  na.samples <- c(na.samples, colnames(eset)[flag])
}

# HERE

if(length(na.samples) != 0) {
  ## cat("Dropping samples that were NA in some feature of interest: ", paste0(na.samples, collapse=","), "\n", sep="")
  anno.df <- anno.df[!(rownames(anno.df) %in% na.samples),]
  eset <- eset[,!(colnames(eset) %in% na.samples)]
}

expr.samples <- expr.samples[!(expr.samples$sanitized %in% na.samples),]
expr.samples <- expr.samples[(expr.samples$sanitized %in% colnames(eset)),]
write.table(file="samples-in-training-and-test-sets.tsv", expr.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

calculate.hr <- function(classifier.res.df) {
  suppressPackageStartupMessages(library("survival"))
  surv.formula <- as.formula(paste("Surv(time, event)", "high.risk", sep=" ~ "))
  cp.fit <- coxph(surv.formula, data = classifier.res.df)  
  sum <- summary(cp.fit)
  hr <- exp(sum$coefficients[1])
  if(hr < 1) { hr <- 1 / hr }
  hr <- format(hr, digits=2)
  hr  
}

plot.km <- function(classifier.res.df, main="") {
  suppressPackageStartupMessages(library("survminer"))
  suppressPackageStartupMessages(library("survival"))
  surv.formula <- as.formula(paste("Surv(time, event)", "high.risk", sep=" ~ "))
  cp.fit <- coxph(surv.formula, data = classifier.res.df)  
  sum <- summary(cp.fit)
  hr <- exp(sum$coefficients[1])
  if(hr < 1) { hr <- 1 / hr }
  hr <- format(hr, digits=2)

  surv.formula <- as.formula(paste("Surv(time, event)", "high.risk", sep=" ~ "))
  km.fit <- survfit(surv.formula, data = classifier.res.df)

  # Don't know why I need to do this.  These values are strings, not objects
  # as they should be, if I don't set explicitly.
  km.fit$call$formula <- surv.formula
  km.fit$call$data <- classifier.res.df
  g <- ggsurvplot(fit=km.fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, main=paste0(main, " HR: ", hr))
}

# Utility function to wrap models
model.entrez <- function(model.name, X, y = NULL, threshold = NULL) {
  supported.models <- c("emc92", "uams17", "uams70", "uams5")
  if(!(model.name %in% supported.models)) {
    cat(paste0(model.name, " is not supported, please specify one of: ", paste(supported.models, collapse=","), "\n"))
    return()
  }
  res <- NULL
  if(model.name == "emc92") {
    res <- emc.92.entrez(X = X, y = y, threshold = threshold)
  } else if(model.name == "uams5") {
    res <- uams.5.entrez(X = X, y = y, threshold = threshold)
  } else if(model.name == "uams17") {
    res <- uams.17.entrez(X = X, y = y, threshold = threshold)
  } else if(model.name == "uams70") {
    res <- uams.70.entrez(X = X, y = y, threshold = threshold)
  }
  return(res)
}

# Break the MMRF data into training and test data sets, stratifying based on risk.
# Set the seed so that our partition is reproducible
set.seed(seed)
# suppressPackageStartupMessages(library(caret))
# training.fraction <- 0.7
# train_ind <- as.vector(createDataPartition(factor(anno.df$PFSbasedHighRisk), p=training.fraction, list = FALSE))

mmrf.flag <- grepl(pattern="MMRF", x=anno.df$ID)

test.eset <- eset[,mmrf.flag]
test.anno <- anno.df[mmrf.flag,]


roc.classifiers <- list()
roc.results <- list()
#cat("Results standardizing by test data\n")
#res <- uams.17.entrez(X = as.matrix(exprs(test.eset)), y = NULL, standardize.using.training.data = FALSE)
#classifier.res <- merge(res$rest$res, test.anno, by="ID")
#ratios.z <- res$ratios

#roc.results[[1]] <- classifier.res
#roc.classifiers[[1]] <- "z-score"

#do.auc(classifier.res, dependent.var="PFSbasedHighRisk")
#do.f1(classifier.res, dependent.var="PFSbasedHighRisk")

#cat("Results standardizing by training data\n")
#res <- uams.17.entrez(X = as.matrix(exprs(test.eset)), y = NULL, standardize.using.training.data = TRUE)
#classifier.res <- merge(res$rest$res, test.anno, by="ID")
#ratios.train <- res$ratios

#roc.results[[2]] <- classifier.res
#roc.classifiers[[2]] <- "train-standard"

#do.auc(classifier.res, dependent.var="PFSbasedHighRisk")
#do.f1(classifier.res, dependent.var="PFSbasedHighRisk")

#png("mmrf-z-score-vs-train-std.png")
#n.classifiers <- length(roc.classifiers)
#palette <- rainbow(n.classifiers, s=0.5)
#aucs <- list()
#rocs <- list()    
#for(i in 1:n.classifiers) {
#  tsv <- roc.results[[i]]
#  response.col <- "PFSbasedHighRisk"
#  predictor.col <- "raw.score"
#  aucs[[i]] <- auc(response=tsv[,response.col], predictor=tsv[,predictor.col])
#  rocs[[i]] <- roc(tsv[,response.col], tsv[,predictor.col])
#}
### Sort ROCs by AUCs
#aucs <- unlist(aucs)
#o <- order(aucs, decreasing=TRUE)
#aucs <- aucs[o]
#rocs <- rocs[o]
#for(i in 1:n.classifiers) {    
#  if(i==1) {
#    plot(rocs[[i]], col=palette[i], main="ROC")
#  } else {
#    plot(rocs[[i]], add=TRUE, col=palette[i])            
#  }
#}

#legends <- unlist(roc.classifiers)
#legends <- legends[o]
#legends <- sapply(1:length(legends), function(i) paste(legends[i], " (AUC: ", format(aucs[[i]], digits=3), ")", sep=""))
#legend(x="bottomright", legend=legends, fill=palette[1:length(legends)])
#d <- dev.off()

#inters <- intersect(names(ratios.z), names(ratios.train))
#ratios.z <- ratios.z[inters]
#ratios.train <- ratios.train[inters]

#png("std-ratios-vs-train-mmrf.png")
#plot(ratios.z, ratios.train, main="Std Dev-scaled UAMS-17 Coefficients", xlab="Z-scored", ylab="Training-based")
#lm.obj <- lm(ratios.z ~ ratios.train)
#r2 <- summary(lm.obj)$adj.r.squared
#abline(lm.obj)
#text(x=min(ratios.z)*0.7, y=max(ratios.train)*0.9, labels=paste0("R2 = ", signif(r2,3)))
#d <- dev.off()

train.eset <- eset[,!mmrf.flag]
train.anno <- anno.df[!mmrf.flag,]

if(any(colnames(exprs(train.eset)) != train.anno$ID)) {
  warning("Training data and annotations are not properly aligned\n")
}

if(any(colnames(exprs(test.eset)) != test.anno$ID)) {
  warning("Test data and annotations are not properly aligned\n")
}

## Create KM curves for the EMC-92, UAMS-17, UAMS-70, and UAMS-5 models.
## For each model, run the model against:
## (1) all data, using the published threshold
## (2) the test data, after optimizing the threshold on the training data
## (3) the test data, using the published threshold
models <- c("emc92", "uams17", "uams70", "uams5")
results <- list()
classifiers <- list()
nxt <- 0
for(i in 1:length(models)) {
  model <- models[i]
  cat(paste0(model, "\n"))

  # Now run the model with training to determine an optimal threshold.
  # Train the classifier on the training data and evaluate on the test data.
  set.seed(seed)
  y <- train.anno$PFSbasedHighRisk
  names(y) <- train.anno$ID
  # Train the classifier -- don't really train; just use the default threshold.  We only 
  # care about AUC/continuous scores for the time being.
  res <- model.entrez(model, X=as.matrix(exprs(train.eset)), y = NULL, threshold = NULL)$obj
  
  # Run the classifier
  trained.res <- model.entrez(model, X=as.matrix(exprs(test.eset)), threshold = res$threshold)
  classifier.res <- merge(trained.res[["res"]], test.anno, by="ID")
  out.file <- paste0(model, "-trained-test.png")
  png(out.file) 
  g <- plot.km(classifier.res, main=paste0(toupper(model)))
  print(g)
  d <- dev.off()

  nxt <- nxt + 1
  results[[nxt]] <- classifier.res
  classifiers[[nxt]] <- paste0(toupper(model))
}

do.auc <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.auc(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$raw.score), digits=3)
}

do.bac <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.bac(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$high.risk), digits=3)
}

do.mcc <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.mcc(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$high.risk), digits=3)
}

do.f1 <- function(data.tsv, dependent.var = "PFSbasedHighRisk") {
  format(calculate.f1(actual=as.logical(data.tsv[,dependent.var]), predicted=data.tsv$high.risk), digits=3)
}

do.hr <- function(data.tsv) {
  format(calculate.hr(data.tsv), digits=2)
}

source("../../novelClassifiers/lm.R")
source("../../novelClassifiers/logit.R")
source("../../novelClassifiers/rf.R")
source("../../novelClassifiers/srf.R")
source("../../performance_metrics/metrics.R")
source("../../novelClassifiers/ensembles.R")
source("../../publishedClassifiers/naive.ensemble/naive.ensemble.R")

rownames(train.anno) <- train.anno$ID
rownames(test.anno) <- test.anno$ID

nxt <- length(results) + 1
set.seed(seed)
res <- naive.ensemble(X = as.matrix(exprs(test.eset)), y = NULL, threshold = NULL)
classifier.res <- merge(res[["res"]], test.anno, by="ID")
model <- "ensemble-geomean"

results[[nxt]] <- classifier.res
classifiers[[nxt]] <- toupper(model)

nxt <- length(results) + 1
set.seed(seed)
res <- ensemble.mean(X = as.matrix(exprs(test.eset)), y = NULL)
classifier.res <- merge(res[["res"]], test.anno, by="ID")
model <- "ensemble-mean"

results[[nxt]] <- classifier.res
classifiers[[nxt]] <- toupper(model)

pfs.results <- results
os.results <- results

pfs.classifiers <- classifiers
os.classifiers <- classifiers

# Run random forest on ISS, sex, age, and expr
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = train.anno[,c("Sex", "ISSstage", "Age")], genomic = NULL)
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = test.anno[,c("Sex", "ISSstage", "Age")], genomic = NULL)

train.mat <- assemble.predictor.matrix(train.data, already.log2.transformed = TRUE)
test.mat <- assemble.predictor.matrix(test.data, already.log2.transformed = TRUE)

for(response in c("PFSbasedHighRisk", "OSbasedHighRisk")) {
  y <- train.anno[,response]
  names(y) <- train.anno$ID
  set.seed(seed)
  rf <- train.random.forest(X = train.mat, y = y)
  res <- test.random.forest(rf, X = test.mat)
  rf.res <- merge(as.data.frame(res), test.anno, by="ID")

  rf.res$high.risk <- rf.res$high.risk == 1

  if(grepl(x=response, pattern="PFS")) {
    nxt <- length(pfs.classifiers) + 1
    pfs.classifiers[[nxt]] <- "RF (ISS+sex+age+expr)"
    pfs.results[[nxt]] <- rf.res
  } else {
    nxt <- length(os.classifiers) + 1
    os.classifiers[[nxt]] <- "RF (ISS+sex+age+expr)"
    os.results[[nxt]] <- rf.res
  }
}

# Run random forest on sex, age, and expr (no ISS)
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = train.anno[,c("Sex", "Age")], genomic = NULL)
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = test.anno[,c("Sex", "Age")], genomic = NULL)

train.mat <- assemble.predictor.matrix(train.data, already.log2.transformed = TRUE)
test.mat <- assemble.predictor.matrix(test.data, already.log2.transformed = TRUE)

for(response in c("PFSbasedHighRisk", "OSbasedHighRisk")) {
  y <- train.anno[,response]
  names(y) <- train.anno$ID
  set.seed(seed)
  rf <- train.random.forest(X = train.mat, y = y)
  res <- test.random.forest(rf, X = test.mat)
  rf.res <- merge(as.data.frame(res), test.anno, by="ID")
  
  rf.res$high.risk <- rf.res$high.risk == 1
  
  if(grepl(x=response, pattern="PFS")) {
    nxt <- length(pfs.classifiers) + 1
    pfs.classifiers[[nxt]] <- "RF (sex+age+expr)"
    pfs.results[[nxt]] <- rf.res
  } else {
    nxt <- length(os.classifiers) + 1
    os.classifiers[[nxt]] <- "RF (sex+age+expr)"
    os.results[[nxt]] <- rf.res
  }
}

# Run random forest on expr (only)
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = NULL, genomic = NULL)
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = NULL, genomic = NULL)

train.mat <- assemble.predictor.matrix(train.data, already.log2.transformed = TRUE)
test.mat <- assemble.predictor.matrix(test.data, already.log2.transformed = TRUE)

for(response in c("PFSbasedHighRisk", "OSbasedHighRisk")) {
  y <- train.anno[,response]
  names(y) <- train.anno$ID
  set.seed(seed)
  rf <- train.random.forest(X = train.mat, y = y)
  res <- test.random.forest(rf, X = test.mat)
  rf.res <- merge(as.data.frame(res), test.anno, by="ID")
  
  rf.res$high.risk <- rf.res$high.risk == 1
  
  if(grepl(x=response, pattern="PFS")) {
    nxt <- length(pfs.classifiers) + 1
    pfs.classifiers[[nxt]] <- "RF (expr)"
    pfs.results[[nxt]] <- rf.res
  } else {
    nxt <- length(os.classifiers) + 1
    os.classifiers[[nxt]] <- "RF (expr)"
    os.results[[nxt]] <- rf.res
  }
}

output.results <- function(classifiers, results, dependent.var = "PFSbasedHighRisk") {
  df.full <- data.frame(Classifier=unlist(classifiers),
                        AUC=unlist(lapply(results, do.auc, dependent.var = dependent.var)))
  
  df.full <- df.full[order(df.full$AUC, decreasing=TRUE),]
  
  write.table(df.full, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  suppressPackageStartupMessages(library(knitr))
  kable(df.full, row.names=FALSE)
  
  library(pROC)
  
  return(list(df.full=df.full))
}

out.pfs <- output.results(pfs.classifiers, pfs.results, dependent.var = "PFSbasedHighRisk")
out.os <- output.results(os.classifiers, os.results, dependent.var = "OSbasedHighRisk")

cat("PFS-based high risk\n")
kable(out.pfs$df.full, row.names=FALSE)

cat("OS-based high risk\n")
kable(out.os$df.full, row.names=FALSE)


df <- matrix(data=0, nrow=length(classifiers), ncol=length(classifiers))
rownames(df) <- unlist(classifiers)
colnames(df) <- unlist(classifiers)


for(i in 1:length(classifiers)) {
  for(j in 1:length(classifiers)) {
    df[i,j] <- cor(results[[i]]$raw.score, results[[j]]$raw.score)
  }
}

write.table(file="correlation-of-classifiers", df, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)