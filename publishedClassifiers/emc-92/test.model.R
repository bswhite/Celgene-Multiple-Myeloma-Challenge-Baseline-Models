source("../common/classifier-base.R")
source("emc-92.R")

use.synapse <- TRUE

if(use.synapse) {
  library(synapseClient)
  synapseLogin()
}

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

cat(paste0("Origin matrix had ", ncol(expr.data), " cols\n"))
expr.data <- expr.data[,!duplicated(cols, fromLast=FALSE)]
cat(paste0("After dropping duplicate samples, matrix has ", ncol(expr.data), " cols\n"))
        
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_BM", replacement="")
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_PB", replacement="")

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
anno.df <- anno.df[,c("SampleID", "PFStimeMonths", "Progression", "PFSbasedHighRisk")]
colnames(anno.df) <- c("ID", "time", "event", "PFSbasedHighRisk")

# Exclude any samples that do have a PFS-based high risk
cat(paste0("Excluding ", length(which(is.na(anno.df$PFSbasedHighRisk))), " of ", nrow(anno.df), " samples without low/high risk\n"))
anno.df <- anno.df[!is.na(anno.df$PFSbasedHighRisk),]
  
inters <- intersect(anno.df$ID, colnames(expr.data))
cat(paste0(length(inters), " columns are shared between the ", nrow(anno.df), " clinically annotated samples and the ", ncol(expr.data), " expression samples\n"))

# Subset the data and annotations to include only samples appearing in both.
anno.df <- anno.df[anno.df$ID %in% inters,]
expr.data <- expr.data[,inters]

# Normalize the data using edgeR
suppressPackageStartupMessages(library(edgeR))

y <- DGEList(counts=expr.data)

# Let's not filter any genes--since the classifier will explicitly pick out those
# of interest.
## Filter non-expressed genes:
## A <- aveLogCPM(y)
## y <- y[A>1,]

# Then normalize and compute log2 counts-per-million with an offset:
# Gordyn uses prior.count = 5 in the above link--I use the default of 0.5
y <- calcNormFactors(y)
expr.tbl <- cpm(y, log=TRUE, prior.count=0.5)
anno.df$ID <- as.character(anno.df$ID)
expr.tbl <- expr.tbl[,anno.df$ID]
eset <- ExpressionSet(as.matrix(expr.tbl))

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

pdf("emc92-kms.pdf")

# Don't specify threshold or training data; this will use the default, published threshold
train.eset <- NULL
train.anno <- NULL
emc.92.published.res <- emc.92.ensg(eset, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = NULL)
classifier.res <- merge(emc.92.published.res[["res"]], anno.df, by="ID")
g <- plot.km(classifier.res, main="Published")
print(g)

# Now run EMC 92 with training to determine an optimal threshold.

# Break the MMRF data into training and test data sets, stratifying based on risk.
# Set the seed so that our partition is reproducible
set.seed(1234)
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

# Train the EMC92 classifier on the training data and evaluate on the test data.
emc.92.trained.res <- emc.92.ensg(test.eset, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = NULL)
classifier.res <- merge(emc.92.trained.res[["res"]], test.anno, by="ID")
g <- plot.km(classifier.res, main="Test Data (Trained Threshold)")
print(g)

# Finally, run EMC 92 on the default threshold, but just on the test data, so we can 
# compare with our trained model.

# Evaluate the classifier on the test data
emc.92.published.test.res <- emc.92.ensg(test.eset, already.log2.transformed = TRUE, train.eset = NULL, train.anno = NULL, threshold = NULL)
classifier.res <- merge(emc.92.published.test.res[["res"]], test.anno, by="ID")
g <- plot.km(classifier.res, main="Test Data (Published Threshold)")
print(g)

d <- dev.off()

