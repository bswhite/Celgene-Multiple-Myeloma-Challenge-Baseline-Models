source("../common/classifier-base.R")
source("../emc-92/emc-92.R")
source("../uams-17/uams-17.R")
source("../uams-70/uams-70.R")

use.synapse <- TRUE

seed <- 1234
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
  
  print(c(a0, a1))
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
  print(table(sigVariedGenes))  
  
  d <- dev.off()
  return(names(means)[sigVariedGenes])
}

# Get the non-synonymous variants
ns.variants.data <- NULL
ns.variants.file <- "MMRF_CoMMpass_IA8b_All_Canonical_NS_Variants.txt"
if(!file.exists(ns.variants.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7450323", downloadFile = TRUE, downloadLocation = ".")
    gz.file <- getFileLocation(obj)
    system(paste0("gunzip ", gz.file))
  }
}
ns.variants.data <- read.table(ns.variants.file, sep="\t", header=TRUE, comment.char="")

data <- NULL
# Get the count data
data.file <- "MMRF_CoMMpass_IA8b_E74GTF_Salmon_Gene_Counts.txt"
if(!file.exists(data.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7450336", downloadFile = TRUE, downloadLocation = ".")
    gz.file <- getFileLocation(obj)
    system(paste0("gunzip ", gz.file))
  }
}

data <- read.table(data.file, sep="\t", header=TRUE)
# TPM data
# data <- read.table("MMRF_CoMMpass_IA8b_E74GTF_Salmon_Gene_TPM.txt", sep="\t", header=TRUE)

# Some of the columns seem to be rederived from the same patient, with names like
# MMRF_2087_1_BM and MMRF_2087_2_BM
# or
# MMRF_1992_1_BM and MMRF_1992_1_PB

# Let's take a single column per patient, favoring "*_1_BM" over "*_1_PB" over "*_2_BM"
# Do this by sorting the columns by name and taking the first non-duplicated.
expr.data <- data[,!(colnames(data) %in% "GENE_ID")]
rownames(expr.data) <- data$GENE_ID

expr.data <- expr.data[,order(colnames(expr.data))]
cols <- colnames(expr.data)
cols <- gsub(cols, pattern="_\\d+_BM", replacement="")
cols <- gsub(cols, pattern="_\\d+_PB", replacement="")

## cat(paste0("Origin matrix had ", ncol(expr.data), " cols\n"))
expr.data <- expr.data[,!duplicated(cols, fromLast=FALSE)]
## cat(paste0("After dropping duplicate samples, matrix has ", ncol(expr.data), " cols\n"))
        
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_BM", replacement="")
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_PB", replacement="")

# Limit the variants to those that are called by at least 2 of 3 callers
flag <- apply(ns.variants.data[,c("SEURAT","MUTECT","STRELKA")], 1, function(row) length(which(row=="true")) >= 2)
ns.variants.data <- ns.variants.data[flag,]

# Let's take a single column per patient, favoring "*_1_BM" over "*_1_PB" over "*_2_BM"
# Do this by sorting the columns by name and taking the first non-duplicated.

ns.variants.data$harmonized.sample <- ns.variants.data$Sample
ns.variants.data$harmonized.sample <- gsub(ns.variants.data$harmonized.sample, pattern="_\\d+_BM", replacement="")
ns.variants.data$harmonized.sample <- gsub(ns.variants.data$harmonized.sample, pattern="_\\d+_PB", replacement="")

## Sort mutations so that those from the patient are ordered with _1_BM first, then _2_BM, _1_PB, etc.
## Then just keep those that are not duplicated--those will be those that are _1_BM or are unique to _2_BM, etc.

## Drop entries that are not mapped to an ENSG gene
## ensg.gene.mutations <- unlist(lapply(ns.variants.data$EFF....GENE, function(x) grepl(x=x, pattern="ENSG")))
## ns.variants.data <- ns.variants.data[ensg.gene.mutations,]

## It appears that are occassionally multiple entries for the same gene/mutation (e.g., based on different)
## transcripts.  Subset to the rows that are unique with respect to the fields we care about
mut.cols <- c("Sample", "X.CHROM", "POS", "REF", "ALT", "GEN.1..AR", "EFF....GENE")
dup <- duplicated(ns.variants.data[,mut.cols])
ns.variants.data <- ns.variants.data[!dup,]

o <- order(do.call(order, ns.variants.data[,mut.cols]))
mut.cols.no.smpl <- c("X.CHROM", "POS", "REF", "ALT", "GEN.1..AR", "EFF....GENE")
## dup <- duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=TRUE) | duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=FALSE)
## head(ns.variants.data[dup,])

dup <- duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=FALSE)
ns.variants.data <- ns.variants.data[!dup,]

## Create a row (ensg gene id) x col (sample) matrix, where the entries are the entries are the highest VAF
## for any mutation in that gene in that sample.
ns.variants.df <- ns.variants.data[,c("harmonized.sample", "EFF....GENE", "GEN.1..AR")]
ns.variants.df <- ns.variants.data[order(ns.variants.data$GEN.1..AR, decreasing=TRUE),]
ns.variants.df <- ns.variants.df[!duplicated(ns.variants.df[,c("harmonized.sample", "EFF....GENE")], fromLast = FALSE),]

suppressPackageStartupMessages(library(reshape2))
gene.by.sample.mutations <- acast(data=ns.variants.df, formula=EFF....GENE ~ harmonized.sample, value.var="GEN.1..AR", fill = 0)

# Pull in the annotation data.  
anno.file <- NULL
if(use.synapse) {
  obj <- synGet(id="syn7488088", downloadFile = TRUE, downloadLocation = ".")
  anno.file <- getFileLocation(obj)
} else {
  anno.file <- "MMRF patient-level clinical and cytogenetic inventory Harmonized.csv"
}

anno.df <- read.table(anno.file, sep=",", header=TRUE)

# This must have columns event, time, and ID.  
# Use PFS annotations to define high risk by setting event = Progression; time = PFStimeMonths; ID = SampleID
anno.df <- anno.df[,c("SampleID", "PFStimeMonths", "Progression", "PFSbasedHighRisk", "Sex", "ISSstage")]
colnames(anno.df) <- c("ID", "time", "event", "PFSbasedHighRisk", "Sex", "ISSstage")

# Exclude any samples that do have a PFS-based high risk
## cat(paste0("Excluding ", length(which(is.na(anno.df$PFSbasedHighRisk))), " of ", nrow(anno.df), " samples without low/high risk\n"))
anno.df <- anno.df[!is.na(anno.df$PFSbasedHighRisk),]

orig.anno.samples <- nrow(anno.df)
orig.expr.samples <- ncol(expr.data)
orig.genomic.samples <- ncol(gene.by.sample.mutations)
rownames(anno.df) <- anno.df$ID
inters <- intersect(rownames(anno.df), colnames(expr.data))
## cat(paste0(length(inters), " columns are shared between the ", nrow(anno.df), " clinically annotated samples and the ", ncol(expr.data), " expression samples\n"))

# Subset the data and annotations to include only samples appearing in both.
anno.df <- anno.df[inters,]
expr.data <- expr.data[,inters]

# Further subset to include only those with genomic data
inters <- intersect(colnames(gene.by.sample.mutations), colnames(expr.data))
## cat(paste0(length(inters), " columns are shared between the ", orig.anno.samples, " clinically annotated sample, the ", orig.expr.samples, " expression samples, and the " , orig.genomic.samples, " genomic samples\n"))
anno.df <- anno.df[inters,]
expr.data <- expr.data[, inters]
gene.by.sample.mutations <- gene.by.sample.mutations[, inters]

# Give the genomic (gene) features names that are distinct from the expression (gene) features
rownames(gene.by.sample.mutations) <- unlist(lapply(rownames(gene.by.sample.mutations), function(x) paste0("mut",x)))

# Normalize the data using edgeR
suppressPackageStartupMessages(library(edgeR))


y <- DGEList(counts=expr.data)

## Filter non-expressed genes (unless they are involved in a classifer):
coefficients <- c(emc.92.get.coefficients(), uams.17.get.coefficients(), uams.70.get.coefficients())

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
probesets <- names(coefficients)
bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'ensembl_gene_id'), 
            filters = 'affy_hg_u133_plus_2', 
            values = probesets, 
            mart = ensembl)
names(bm) <- c("PROBE", "GENE")
coefficient.genes <- bm$GENE

A <- aveLogCPM(y)
y <- y[(A>1) | (rownames(y$counts) %in% coefficient.genes),]

# Then normalize and compute log2 counts-per-million with an offset:
# Gordyn uses prior.count = 5 in the above link--I use the default of 0.5
y <- calcNormFactors(y)
expr.tbl <- cpm(y, log=TRUE, prior.count=0.5)
anno.df$ID <- as.character(anno.df$ID)
expr.tbl <- expr.tbl[,anno.df$ID]
eset <- ExpressionSet(as.matrix(expr.tbl))

# Only test those genes that are mutated in some sample
flag <- unlist(apply(gene.by.sample.mutations, 1, function(row) any(row>0)))
genomic.gene.subset <- rownames(gene.by.sample.mutations)[flag]

## Only test genes that are mutated at a rate > 2%
flag <- unlist(apply(gene.by.sample.mutations, 1, function(row) ( length(which(row>0))*1.0/length(row) > 0.02 )))
genomic.gene.subset <- rownames(gene.by.sample.mutations)[flag]
empirically.highly.mutated.genes <- unlist(rownames(gene.by.sample.mutations)[flag], function(x) gsub(x=x, pattern="mut", replacement=""))

translocation.targets <- c("FGFR3", "MMSET", "WHSC1", "CCND3", "CCND1", "MAF", "MAFB")

highly.mutated.genes <- c("CDKN2C", "RB1", "CCND1", "CDKN2A", "NRAS", "KRAS", "BRAF", "MYC", "TRAF3", "CYLD", "DKK1", "FRZB", "DNAH5", "XBP1", "BLIMP1", "PRDM1", "IRF4", "TP53", "MRE11A", "PARP1", "DIS3", "FAM46C", "LRRK2", "KDM6A", "UTX", "MLL", "MMSET", "WHSC1", "HOXA9", "KDM6B", "IGLL5")

sym.to.ensg <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
                     filters = 'hgnc_symbol', 
                     values = c(translocation.targets, highly.mutated.genes, empirically.highly.mutated.genes), 
                     mart = ensembl)
names(sym.to.ensg) <- c("SYMBOL", "GENE")

## When looking at the expression data, check genes that are highly mutated/translocated or are used by one
## of the other classifiers, and are expressed
expr.gene.subset <- intersect(rownames(eset), sym.to.ensg$GENE)
expr.gene.subset <- c(expr.gene.subset, coefficient.genes)
expr.gene.subset <- c(expr.gene.subset, get.variable.genes(expr.data))
expr.gene.subset <- intersect(expr.gene.subset, rownames(eset))
cat(paste0("Number of expression features: ", length(expr.gene.subset), "\n"))

## Restrict to highly variable genes

## Use every expression feature
## expr.gene.subset <- rownames(eset)

## Use every genomic feature
## genomic.gene.subset <- rownames(gene.by.sample.mutations)

## Let's exclude IGK, IGH, IGLV (but not IGGL5), and TTN
exclude.patterns <- c("IGK", "IGH", "IGLV", "TTN", "MUC")
exclude <- rep(FALSE, length(genomic.gene.subset))
for(pattern in exclude.patterns) {
  exclude <- exclude | grepl(x=genomic.gene.subset, pattern=pattern)
}
genomic.gene.subset <- genomic.gene.subset[!exclude]

genomic.gene.subset <- intersect(genomic.gene.subset, rownames(gene.by.sample.mutations))

# Finally, make sure we don't have any NAs in any of our features.  If we do, drop the corresponding sample.
na.samples <- c()
flag <- unlist(apply(anno.df[,c("Sex", "ISSstage")], 1, function(row) any(is.na(row))))
if(any(flag)) {
  na.samples <- c(na.samples, rownames(anno.df)[flag])
}

flag <- unlist(apply(gene.by.sample.mutations[genomic.gene.subset,], 2, function(col) any(is.na(col))))
if(any(flag)) {
  na.samples <- c(na.samples, colnames(gene.by.sample.mutations)[flag])
}

flag <- unlist(apply(exprs(eset)[expr.gene.subset,], 2, function(col) any(is.na(col))))
if(any(flag)) {
  na.samples <- c(na.samples, colnames(eset)[flag])
}

if(length(na.samples) != 0) {
  ## cat("Dropping samples that were NA in some feature of interest: ", paste0(na.samples, collapse=","), "\n", sep="")
  anno.df <- anno.df[!(rownames(anno.df) %in% na.samples),]
  eset <- eset[,!(colnames(eset) %in% na.samples)]
  gene.by.sample.mutations <- gene.by.sample.mutations[,!(colnames(gene.by.sample.mutations) %in% na.samples)]
}

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
model.ensg <- function(model.name, test.eset, test.anno = NULL, already.log2.transformed = FALSE, train.eset = NULL, train.anno = NULL, threshold = NULL) {
  supported.models <- c("emc92", "uams17", "uams70")
  if(!(model.name %in% supported.models)) {
    cat(paste0(model.name, " is not supported, please specify one of: ", paste(supported.models, collapse=","), "\n"))
    return()
  }
  res <- NULL
  if(model.name == "emc92") {
    res <- emc.92.ensg(test.eset = test.eset, test.anno = test.anno, already.log2.transformed = already.log2.transformed, train.eset = train.eset, train.anno = train.anno, threshold = threshold)
  } else if(model.name == "uams17") {
    res <- uams.17.ensg(test.eset = test.eset, test.anno = test.anno, already.log2.transformed = already.log2.transformed, train.eset = train.eset, train.anno = train.anno, threshold = threshold)
  } else if(model.name == "uams70") {
    res <- uams.70.ensg(test.eset = test.eset, test.anno = test.anno, already.log2.transformed = already.log2.transformed, train.eset = train.eset, train.anno = train.anno, threshold = threshold)
  }
  return(res)
}

# Break the MMRF data into training and test data sets, stratifying based on risk.
# Set the seed so that our partition is reproducible
set.seed(seed)
suppressPackageStartupMessages(library(caret))
training.fraction <- 0.7
train_ind <- as.vector(createDataPartition(factor(anno.df$PFSbasedHighRisk), p=training.fraction, list = FALSE))

train.eset <- eset[,train_ind]
train.anno <- anno.df[train_ind,]
test.eset <- eset[,-train_ind]
test.anno <- anno.df[-train_ind,]

if(any(colnames(exprs(train.eset)) != train.anno$ID)) {
  warning("Training data and annotations are not properly aligned\n")
}

if(any(colnames(exprs(test.eset)) != test.anno$ID)) {
  warning("Test data and annotations are not properly aligned\n")
}

## Create KM curves for the EMC-92, UAMS-17, and UAMS-70 models.
## For each model, run the model against:
## (1) all data, using the published threshold
## (2) the test data, after optimizing the threshold on the training data
## (3) the test data, using the published threshold
models <- c("emc92", "uams17", "uams70")
results <- list()
classifiers <- list()
for(i in 1:length(models)) {
  model <- models[i]
  ## cat(paste0(model, "\n"))

  # Don't specify threshold or training data; this will use the default, published threshold
  set.seed(seed)
  published.res <- model.ensg(model, test.eset = eset, test.anno = anno.df, already.log2.transformed = TRUE, train.eset = NULL, train.anno = NULL, threshold = NULL)
  classifier.res <- merge(published.res[["res"]], anno.df, by="ID")
  out.file <- paste0(model, "-published-all.png")
  png(out.file) 
  g <- plot.km(classifier.res, main="Published")
  print(g)
  d <- dev.off()

  # Now run the model with training to determine an optimal threshold.

  # Train the classifier on the training data and evaluate on the test data.
  set.seed(seed)
  trained.res <- model.ensg(model, test.eset = test.eset, test.anno = test.anno, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = NULL)
  classifier.res <- merge(trained.res[["res"]], test.anno, by="ID")
  out.file <- paste0(model, "-trained-test.png")
  png(out.file) 
  g <- plot.km(classifier.res, main=paste0(toupper(model), " (Trained Threshold)"))
  print(g)
  d <- dev.off()

  results[[(2*(i-1))+1]] <- classifier.res
  classifiers[[(2*(i-1))+1]] <- paste0(toupper(model), " (Trained)")
  

  # Finally, run classifier using the default threshold, but just on the test data, so we can 
  # compare with our trained model.

  # Evaluate the classifier on the test data
  set.seed(seed)
  published.test.res <- model.ensg(model, test.eset = test.eset, test.anno = test.anno, already.log2.transformed = TRUE, train.eset = NULL, train.anno = NULL, threshold = NULL)
  classifier.res <- merge(published.test.res[["res"]], test.anno, by="ID")
  out.file <- paste0(model, "-published-test.png")
  png(out.file) 
  g <- plot.km(classifier.res, main=paste0(toupper(model), " (Published Threshold)"))
  print(g)
  d <- dev.off()

  results[[2*i]] <- classifier.res
  classifiers[[2*i]] <- paste0(toupper(model), " (Published)")

}

do.auc <- function(data.tsv) {
  format(calculate.auc(actual=as.logical(data.tsv$PFSbasedHighRisk), predicted=data.tsv$raw.score), digits=3)
}

do.bac <- function(data.tsv) {
  format(calculate.bac(actual=as.logical(data.tsv$PFSbasedHighRisk), predicted=data.tsv$high.risk), digits=3)
}

do.mcc <- function(data.tsv) {
  format(calculate.mcc(actual=as.logical(data.tsv$PFSbasedHighRisk), predicted=data.tsv$high.risk), digits=3)
}

do.f1 <- function(data.tsv) {
  format(calculate.f1(actual=as.logical(data.tsv$PFSbasedHighRisk), predicted=data.tsv$high.risk), digits=3)
}

do.hr <- function(data.tsv) {
  format(calculate.hr(data.tsv), digits=2)
}

source("../../novelClassifiers/rf.R")
source("../../novelClassifiers/srf.R")
source("../../performance_metrics/metrics.R")

rownames(train.anno) <- train.anno$ID
rownames(test.anno) <- test.anno$ID

##y <- DGEList(counts=expr.data)

### Renormalize the data after excluding lowly expressed genes
##A <- aveLogCPM(y)
##y <- y[A>1,]

### Then normalize and compute log2 counts-per-million with an offset:
### Gordyn uses prior.count = 5 in the above link--I use the default of 0.5
##y <- calcNormFactors(y)
##expr.tbl <- cpm(y, log=TRUE, prior.count=0.5)
##expr.tbl <- expr.tbl[,anno.df$ID]
##eset <- ExpressionSet(as.matrix(expr.tbl))

train.eset <- eset[,train_ind]
train.anno <- anno.df[train_ind,]
train.genomic <- gene.by.sample.mutations[,train_ind]

test.eset <- eset[,-train_ind]
test.anno <- anno.df[-train_ind,]
test.genomic <- gene.by.sample.mutations[,-train_ind]

train.genomic[train.genomic > 0] <- 1
test.genomic[test.genomic > 0] <- 1

## Run RF with ISS+mut+expr
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = train.anno[,c("Sex", "ISSstage")], genomic = train.genomic[genomic.gene.subset,])
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = test.anno[,c("Sex", "ISSstage")], genomic = test.genomic[genomic.gene.subset,])
set.seed(seed)
rf <- train.random.forest(data = train.data, labels = train.anno, already.log2.transformed = TRUE)

res <- run.random.forest(rf, data = test.data, labels = test.anno, already.log2.transformed = TRUE)
rf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(rf.res$high.risk)
rf.res$high.risk <- r == 1

classifier.name <- "RF (ISS+mut+expr)"
classifiers[[7]] <- classifier.name
results[[7]] <- rf.res

out.file <- "RF-ISS-mut-expr-trained-test.png"
png(out.file) 
g <- plot.km(rf.res, main=classifier.name)
print(g)
d <- dev.off()

train.surv <- train.anno[,c("time", "event")]
colnames(train.surv) <- c("time", "status")
test.surv <- test.anno[,c("time", "event")]
colnames(test.surv) <- c("time", "status")

test.labels = test.anno$PFSbasedHighRisk
names(test.labels) <- rownames(test.anno)

train.labels = train.anno$PFSbasedHighRisk
names(train.labels) <- rownames(train.anno)

set.seed(seed)
srf <- train.surv.random.forest(data = train.data, surv.anno = train.surv, labels = train.labels, already.log2.transformed = TRUE)
res <- run.surv.random.forest(srf, data = test.data, surv.anno = test.surv, labels = test.labels, already.log2.transformed = TRUE)
srf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(srf.res$high.risk)
srf.res$high.risk <- r == 1

classifier.name <- "SRF (ISS+mut+expr)"
classifiers[[8]] <- classifier.name
results[[8]] <- srf.res

out.file <- "SRF-ISS-mut-expr-trained-test.png"
png(out.file) 
g <- plot.km(srf.res, main=classifier.name)
print(g)
d <- dev.off()

## Run RF with ISS+mut
train.data <- list(eset = NULL, clinical = train.anno[,c("Sex", "ISSstage")], genomic = train.genomic[genomic.gene.subset,])
test.data <- list(eset = NULL, clinical = test.anno[,c("Sex", "ISSstage")], genomic = test.genomic[genomic.gene.subset,])
set.seed(seed)
rf <- train.random.forest(data = train.data, labels = train.anno, already.log2.transformed = TRUE)

res <- run.random.forest(rf, data = test.data, labels = test.anno, already.log2.transformed = TRUE)
rf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(rf.res$high.risk)
rf.res$high.risk <- r == 1

classifier.name <- "RF (ISS+mut)"
classifiers[[9]] <- classifier.name
results[[9]] <- rf.res

out.file <- "RF-ISS-mut-trained-test.png"
png(out.file) 
g <- plot.km(rf.res, main=classifier.name)
print(g)
d <- dev.off()

set.seed(seed)
srf <- train.surv.random.forest(data = train.data, surv.anno = train.surv, labels = train.labels, already.log2.transformed = TRUE)
res <- run.surv.random.forest(srf, data = test.data, surv.anno = test.surv, labels = test.labels, already.log2.transformed = TRUE)
srf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(srf.res$high.risk)
srf.res$high.risk <- r == 1

classifier.name <- "SRF (ISS+mut)"
classifiers[[10]] <- classifier.name
results[[10]] <- srf.res

out.file <- "SRF-ISS-mut-trained-test.png"
png(out.file) 
g <- plot.km(srf.res, main=classifier.name)
print(g)
d <- dev.off()


## Run RF with ISS+expr
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = train.anno[,c("Sex", "ISSstage")], genomic = NULL)
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = test.anno[,c("Sex", "ISSstage")], genomic = NULL)
set.seed(seed)
rf <- train.random.forest(data = train.data, labels = train.anno, already.log2.transformed = TRUE)

res <- run.random.forest(rf, data = test.data, labels = test.anno, already.log2.transformed = TRUE)
rf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(rf.res$high.risk)
rf.res$high.risk <- r == 1

classifiers[[11]] <- "RF (ISS+expr)"
results[[11]] <- rf.res

out.file <- "RF-ISS-expr-trained-test.png"
png(out.file) 
g <- plot.km(rf.res, main="RF")
print(g)
d <- dev.off()

set.seed(seed)
srf <- train.surv.random.forest(data = train.data, surv.anno = train.surv, labels = train.labels, already.log2.transformed = TRUE)
res <- run.surv.random.forest(srf, data = test.data, surv.anno = test.surv, labels = test.labels, already.log2.transformed = TRUE)
srf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(srf.res$high.risk)
srf.res$high.risk <- r == 1

classifier.name <- "SRF (ISS+expr)"
classifiers[[12]] <- classifier.name
results[[12]] <- srf.res

out.file <- "SRF-iss-expr-trained-test.png"
png(out.file) 
g <- plot.km(srf.res, main=classifier.name)
print(g)
d <- dev.off()


## Run RF with expr
train.data <- list(eset = train.eset[expr.gene.subset,], clinical = NULL, genomic = NULL)
test.data <- list(eset = test.eset[expr.gene.subset,], clinical = NULL, genomic = NULL)
set.seed(seed)
rf <- train.random.forest(data = train.data, labels = train.anno, already.log2.transformed = TRUE)

res <- run.random.forest(rf, data = test.data, labels = test.anno, already.log2.transformed = TRUE)
rf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(rf.res$high.risk)
rf.res$high.risk <- r == 1

classifier.name <- "RF (expr)"
classifiers[[13]] <- classifier.name
results[[13]] <- rf.res

out.file <- "RF-expr-trained-test.png"
png(out.file) 
g <- plot.km(rf.res, main=classifier.name)
print(g)
d <- dev.off()

set.seed(seed)
srf <- train.surv.random.forest(data = train.data, surv.anno = train.surv, labels = train.labels, already.log2.transformed = TRUE)
res <- run.surv.random.forest(srf, data = test.data, surv.anno = test.surv, labels = test.labels, already.log2.transformed = TRUE)
srf.res <- merge(as.data.frame(res), test.anno, by="ID")

r <- as.numeric(srf.res$high.risk)
srf.res$high.risk <- r == 1

classifier.name <- "SRF (expr)"
classifiers[[14]] <- classifier.name
results[[14]] <- srf.res

out.file <- "SRF-expr-trained-test.png"
png(out.file) 
g <- plot.km(srf.res, main=classifier.name)
print(g)
d <- dev.off()


df <- data.frame(Classifier=unlist(classifiers),
                  AUC=unlist(lapply(results, do.auc)),
                  BAC=unlist(lapply(results, do.bac)),
                  F1=unlist(lapply(results, do.f1)),
                  MCC=unlist(lapply(results, do.mcc)),
                  HR=unlist(lapply(results, do.hr)))

df <- df[order(df$AUC, decreasing=TRUE),]

write.table(df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

suppressPackageStartupMessages(library(knitr))
kable(df, row.names=FALSE)


library(pROC)

# "Trained" and "Published" models give the same AUC--so just include one of them
flag <- !grepl(x=classifiers, pattern="Trained")
roc.classifiers <- unlist(lapply(classifiers[flag], function(x) gsub(x=x, pattern=" (Published)", replacement="", fixed=TRUE)))
roc.results <- results[flag]

# Also, exclude SRF
flag <- !grepl(x=roc.classifiers, pattern="SRF")
roc.classifiers <- roc.classifiers[flag]
roc.results <- roc.results[flag]

# Exclude some others just for space
## flag <- grepl(x=roc.classifiers, pattern="(expr)", fixed=TRUE) | grepl(x=roc.classifiers, pattern="(ISS+mut)", fixed=TRUE)
flag <- grepl(x=roc.classifiers, pattern="(ISS+mut)", fixed=TRUE)
roc.classifiers <- roc.classifiers[!flag]
roc.results <- roc.results[!flag]


df <- data.frame(Classifier=unlist(roc.classifiers),
                 AUC=unlist(lapply(roc.results, do.auc)),
                 BAC=unlist(lapply(roc.results, do.bac)),
                 F1=unlist(lapply(roc.results, do.f1)),
                 MCC=unlist(lapply(roc.results, do.mcc)),
                 HR=unlist(lapply(roc.results, do.hr)))

df <- df[order(df$AUC, decreasing=TRUE),]

kable(df, row.names=FALSE)
write.table(df, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

n.classifiers <- length(roc.classifiers)
palette <- rainbow(n.classifiers, s=0.5)
aucs <- list()
rocs <- list()    
for(i in 1:n.classifiers) {
  tsv <- roc.results[[i]]
  response.col <- "PFSbasedHighRisk"
  predictor.col <- "raw.score"
  aucs[[i]] <- auc(response=tsv[,response.col], predictor=tsv[,predictor.col])
  rocs[[i]] <- roc(tsv[,response.col], tsv[,predictor.col])
}
## Sort ROCs by AUCs
aucs <- unlist(aucs)
o <- order(aucs, decreasing=TRUE)
aucs <- aucs[o]
rocs <- rocs[o]
for(i in 1:n.classifiers) {    
  if(i==1) {
    plot(rocs[[i]], col=palette[i], main="ROC")
  } else {
    plot(rocs[[i]], add=TRUE, col=palette[i])            
  }
}

legends <- unlist(roc.classifiers)
legends <- legends[o]
legends <- sapply(1:length(legends), function(i) paste(legends[i], " (AUC: ", format(aucs[[i]], digits=3), ")", sep=""))
legend(x="bottomright", legend=legends, fill=palette[1:length(legends)])

# TODO
# - fix implementation of uams70--the 1/-1 one--to average
# - check in uams70
# - uams 5
# - check in uams5

# - update git hub
# - call iss
# - call emc92 + iss

# - rerun
# - check in everything

# - upload table with HR and AUC
# - upload ROC
# - upload KM curves

# - make into slides