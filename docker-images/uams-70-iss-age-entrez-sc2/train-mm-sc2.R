#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("maxstat"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("mice"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("maxstat"))

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

synapseLogin()

## docker build -t docker.synapse.org/syn7237116/uams-70-entrez-sc2:version2 .
## docker images
## docker login docker.synapse.org
## docker push docker.synapse.org/syn7237116/uams-70-entrez-sc2:version2

base.dir <- "./"

download.dir <- paste0(base.dir, "download/")
dir.create(download.dir)

## Read in the annotation file, which indicates the files
## holding the training data.
chal_data_table <- synTableQuery('select id,name from syn9763946')
chal_data_df <- chal_data_table@values
file <- "sc2_Training_ClinAnnotations.csv"
chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern=file))
sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=download.dir)})
training.validation.file <- paste0(download.dir, file)
anno.tbl <- read.table(training.validation.file, sep=",", header=TRUE)

## Read in the clinical annotation file that lists the samples and what data
## are available for each in the columns:
data.set.cols <- colnames(anno.tbl)[grepl(pattern="File", colnames(anno.tbl)) & !grepl(pattern="Sampl", colnames(anno.tbl))]

## The mapping between the Patient/row of the clinical annotation file and
## the identifier in each data file is provided by the corresponding ID
## columns of the clinical annotation file.
## These have the same names as the data set columns, but with "SamplId"
## appended.
data.set.patient.id.cols <- unlist(lapply(data.set.cols,
                                          function(str) paste0(str, "SamplId")))

## Restrict to gene-level expression data sets
expression.data.set.flag <- grepl(pattern="geneLevelExp", data.set.cols)
data.set.cols <- data.set.cols[expression.data.set.flag]
data.set.patient.id.cols <- data.set.patient.id.cols[expression.data.set.flag]

## Only download the training data that is pertinent to this subchallenge.
## Subset to the files we are interested in.
all.data.sets <- unique(na.omit(as.vector(as.matrix(anno.tbl[,data.set.cols]))))
chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
## Download the data and put it in the base directory of the docker image,
## which is where it would be mounted in the challenge.
sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=download.dir)})

## Read in the model implementation.
source(paste0(base.dir, "R/classifier-base.R"))
source(paste0(base.dir, "R/uams-70.R"))

## Read in each of the data sets and create:
## 1) A list holding each of the data sets.
## 2) A list with an entry for each data set holding a mapping from the
##    sample name in the clinical annotation column and the sample name
##    in the expression data set.

## data.sets will eventually hold only those genes that (1) are common
## across all data sets _and_ (2) are in the model.
## Whereas, data.sets.all.common.genes will hold all genes that are
## common across all data sets.
data.sets <- list()
sample.name.mappings <- list()

for(col.indx in 1:length(data.set.cols)) {
    for(data.set in na.omit(unique(anno.tbl[,data.set.cols[col.indx]]))) {
        cat(paste0("Reading in data.set ", data.set, "\n"))
        file <- paste0(download.dir, data.set)
        tbl <- as.data.frame(fread(file, header=TRUE))

        ## The first column is the gene identifier--make it the row
        ## name and drop it from the table
        rownames(tbl) <- tbl[,1]
        tbl <- tbl[,-1]
        data.sets[[data.set]] <- tbl

        ## Extract the sample name mappings from the annotation file
        patient.id.col <- data.set.patient.id.cols[col.indx]
        flag <- !is.na(anno.tbl[,data.set.cols[col.indx]]) & (anno.tbl[,data.set.cols[col.indx]] == data.set)
        map <- anno.tbl[flag,c("Patient", patient.id.col)]
        colnames(map) <- c("patient.id", "sample.id")
        
        sample.name.mappings[[data.set]] <- map
    }
}

## Combine multiple sample columns into a single patient column by
## taking the average of the sample columns
## Columns (i.e., samples) are assumed to have names in sample.to.patient.map$from.
## They will be renamed as sample.to.patient.map$to.
combine_samples_2_patient <- function(expr, sample.to.patient.map) {
    map <- subset(sample.to.patient.map, from %in% colnames(expr))
    rownames(map) <- map$from
    expr <- expr[, colnames(expr) %in% map$from]
    map <- map[colnames(expr),]
    expr <- as.matrix(t(aggregate.Matrix(t(expr), groupings=list(map$to), fun = "mean")))
    expr
}

## Patients listed in the annotation file may have 0, 1, or >1 corresponding
## samples in the data file.  If there are multiple, take the mean across
## samples (NB: since most of these are log expression values, this corresponds
## to the geometric mean in real space).
## Also, subset both the mapping tables and the data to have the
## intersection of patients in both.
for(data.set in names(sample.name.mappings)) {
    tbl <- data.sets[[data.set]]
    map <- sample.name.mappings[[data.set]]

    map <- na.omit(map)

    ## If the mapping from patient to sample id is 1-to-many, split it into multiple
    ## rows, each of which maps the patient to a single sample id.
    map <- ldply(1:nrow(map),
                 .fun = function(i) {
                     sample.id <- unlist(strsplit(as.character(map$sample.id[i]), split=";[ ]*"))
                     df <- data.frame(patient.id = map$patient.id[i], sample.id = sample.id)
                     df
                 })
    
    map <- map[,c("patient.id", "sample.id")]
    colnames(map) <- c("to", "from")
    tbl <- combine_samples_2_patient(tbl, map)

    ## At this point, we have renamed and combined the sample ids columns to patient ids
    ## so that the map is just the identity
    map <- data.frame(patient.id = colnames(tbl), sample.id = colnames(tbl))
    
    data.sets[[data.set]] <- tbl
    sample.name.mappings[[data.set]] <- map
}


## Assemble the survival data.
map <- do.call("rbind", sample.name.mappings)
rownames(map) <- NULL
map <- unique(map)
if(any(duplicated(map$patient.id))) {
    warning("Was not expecting any patient IDs to be duplicated\n")
}

## The sample IDs are scattered across multiple columns, put them
## in the single "sample.id" column
anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)

## Subset the expression data sets to only have the genes common across
## all data sets.

all.genes <- lapply(data.sets, rownames)
common.genes <- Reduce(intersect, all.genes)

## Also, intersect with genes present in the validation data (which are stored on synapse in file overlappingEntrezIds).
## This will ensure that we train our model on genes that are guaranteed to be present in the all data sets (training
## and validation).
synId <- "syn10792130"
obj <- synGet(synId, downloadLocation=download.dir)
overlapping.genes <- read.table(getFileLocation(obj), sep="\t", header=TRUE)

common.genes <- intersect(common.genes, overlapping.genes$entrezID)

## Further subset the expression data to only include those genes in
## the UAMS-70 model.

## Get the Entrez genes in the model.  This is defined in R/uams-70.R
model.entrezs <- uams.70.probe.to.entrez.mapping()
model.entrezs <- na.omit(model.entrezs$INDEX)
common.model.genes <- intersect(common.genes, model.entrezs)

data.sets.all.common.genes <- data.sets

## Subset all of the data sets
for(data.set in names(data.sets)) {
    data.sets.all.common.genes[[data.set]] <- data.sets.all.common.genes[[data.set]][as.character(common.genes),]
    data.sets[[data.set]] <- data.sets[[data.set]][as.character(common.model.genes),]
}

## Z-score each gene, separately for each data set.
## NB: scale normalizes the _columns_ of the matrix
## NB: also we are not scaling the genes in data.sets.all.common.genes,
##     which will be used for imputation.
for(data.set in names(data.sets)) {
    data.sets[[data.set]] <- t(scale(t(data.sets[[data.set]]), center = TRUE, scale = TRUE))
}

## Combine all of the expression matrices into a single matrix and
## ensure they are ordered consistently with the survival data.
expr <- do.call("cbind", data.sets)
expr <- expr[, anno.tbl$sample.id]

expr.all.common.genes <- do.call("cbind", data.sets.all.common.genes)
expr.all.common.genes <- expr.all.common.genes[, anno.tbl$sample.id]

## The ISS column in the clinical annotation file--it has entries,
## 1, 2, 3, or NA.
iss.col <- "D_ISS"

## Train the UAMS-70 model.
y <- anno.tbl[,c("sample.id", "D_PFS", "D_PFS_FLAG", iss.col)]

print(table(is.na(y$D_ISS)))

## FALSE  TRUE 
##  1567    53 

all(anno.tbl$Patient %in% colnames(expr))

## NB: when imputing need to use all genes, not just those in UAMS model.
## Also, don't zscore.  Hence, use the genes in expr.all.common.genes.

imp.cols <- c("Study", "D_Age", "D_Gender", iss.col)
imp.matrix <- as.data.frame(t(expr.all.common.genes))
## ## Just use this for testing expediency:
## imp.matrix <- as.data.frame(t(expr))
imp.matrix$Patient <- rownames(imp.matrix)
imp.matrix <- merge(imp.matrix, anno.tbl[, c("Patient", imp.cols)], by = "Patient")
rownames(imp.matrix) <- imp.matrix$Patient
imp.matrix <- imp.matrix[, colnames(imp.matrix) != "Patient"]
imp.matrix[, iss.col] <- factor(imp.matrix[, iss.col])

cat("Imputing via mice\n")
tempData <- mice(imp.matrix,m=1,maxit=5,meth='pmm',seed=500)
cat("Done imputing via mice\n")
comp <- complete(tempData)

colnames(y) <- c("ID", "time", "event", "ISS")
rownames(y) <- y$ID
y <- na.omit(y)

## Merge in the imputed ISS and age
y <- y[, !(colnames(y) == "ISS")]
comp$Patient <- rownames(comp)
tmp <- comp[, c("Patient", "D_Age", "D_Gender", "D_ISS")]
colnames(tmp) <- c("Patient", "age", "gender", "iss")
y <- merge(y, tmp, by.x = "ID", by.y = "Patient")

rownames(y) <- y$ID
common <- intersect(rownames(y), colnames(expr))
y <- y[common, ]
expr <- expr[, common]

## Invoke UAMS70 model to get scores.  We will recompute threshold after incorporating
## age, gender, and ISS.  So, just use a dummy threshold to avoid optimizing it.
res <- uams.70.entrez(X = as.matrix(expr), y = y, threshold = 0)

surv.features <- merge(res$res[, c("ID", "raw.score")], y[, c("ID", "age", "gender", "iss", "time", "event")], by = "ID")

surv.formula <- as.formula("Surv(time, event) ~ raw.score + age + gender + iss")
fit <- coxph(surv.formula, data = surv.features)

## Calculate the risk
risk <- predict(fit, newdata = surv.features, type="risk")

tmp <- surv.features
tmp$risk <- risk
tmp$high.risk <- unlist(apply(tmp[, c("event", "time")], 1, function(row) ifelse((row[1] == 1) && (row[2] < 18*(365.25/12)), "high", "low")))
##library(ggplot2)
##g <- ggplot(data = tmp, aes(x = high.risk, y = risk))
##g <- g + geom_boxplot()

maxstat.formula <- as.formula("Surv(time, event) ~ risk")
mt <- maxstat.test(maxstat.formula, tmp, smethod="LogRank", pmethod="none")
threshold <- unname(mt$estimate)

tmp$predicted.high.risk <- unlist(lapply(tmp$risk > threshold, function(x) ifelse(x, "high", "low")))
table(tmp$predicted.high.risk, tmp$high.risk)

## Save the trained state of the model, which consistts of
## (1) "threshold" -- the threshold for discretizing the continuous prediction into high vs low risk.
## (2) "mapping" -- a mapping from probe ids to entrez ids.  We will need this mapping during the
## validation phase where we don't have access to the internet.

mapping <- uams.70.probe.to.entrez.mapping()
model.state <- list("threshold" = threshold, "mapping" = mapping, "cox.fit" = fit)

model.state.metadata.file <- paste0(base.dir, "/", "model-state-metadata.Rd")
save(model.state, file=model.state.metadata.file)

cat("Successfully trained model\n")
