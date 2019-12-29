#!/usr/bin/env Rscript

## R library dependencies.  These must be installed in the Docker image.
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("randomForestSRC"))

print(sessionInfo())

use.sequestered.data <- FALSE

if(use.sequestered.data) {
  suppressPackageStartupMessages(library("synapseClient"))
  synapseLogin()
}

## This R script assumes that it runs in the root directory of the Docker image and that
## the test data are mounted in ./test-data,
## the output should be written to ./output,
## and that the entire directory structure of the submitted Docker image
## (e.g., an R object encapsulating trained modeler state) is mounted at ./
docker.image.base <- "./"
test.dir <- "./test-data/"
output.dir <- "./output/"

## Read in the trained model state.
model.state.metadata.file <- paste0(docker.image.base, "model-state-metadata.Rd")
load(model.state.metadata.file)
## model.state is assumed defined in model.state.metadata.file.

## It is assumed to be a list with three entries:
## "model"--the trained random forest
## "target.gene.entrez"--the list of entrez genes that are features (along with D_ISS, D_Age, and D_Gender)
## in the model.
## "y.genes"--entrez ids of genes on the y chromosome
trained.rf <- model.state[["model"]]
target.gene.entrez <- model.state[["target.gene.entrez"]]
y.genes <- model.state[["y.genes"]]
cat(paste0("Class of trained.rf = ", class(trained.rf), "\n"))

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

read.annotation.file <- function() {

    ## This is a kludge to run this without Docker.  Download only the files we are
    ## interested in from the sequester site (at synId = syn9763945)
    chal_data_df <- NULL
    
    if(use.sequestered.data) {
        chal_data_table <- synTableQuery('select id,name from syn9763945')
        chal_data_df <- chal_data_table@values
        
        ## Download the annotation file
        chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc2_Validation_ClinAnnotations.csv"))
        sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=test.dir)})
    }

    ## Read in the annotation file
    validation.annotation.file <- paste0(test.dir, "sc2_Validation_ClinAnnotations.csv")
    anno.tbl <- read.table(validation.annotation.file, sep=",", header=TRUE)
    anno.tbl
}

process.expression.data.sets <- function(anno.tbl) {

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
    
    if(use.sequestered.data) {
        chal_data_table <- synTableQuery('select id,name from syn9763945')
        chal_data_df <- chal_data_table@values
        chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
        ## Download the data and put it in the base directory of the docker image,
        ## which is where it would be mounted in the challenge.
        sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=test.dir)})
    }
    
    ## Read in each of the data sets and create:
    ## 1) A list holding each of the data sets.
    ## 2) A list with an entry for each data set holding a mapping from the
    ##    sample name in the clinical annotation column and the sample name
    ##    in the expression data set.
    data.sets <- list()
    sample.name.mappings <- list()

    ## Iterate over each column holding data
    for(col.indx in 1:length(data.set.cols)) {

        ## Iterate over each data set with data in that column
        for(data.set in na.omit(unique(anno.tbl[,data.set.cols[col.indx]]))) {

            ## Read in the data for this data set
            file <- paste0(test.dir, data.set)
            cat(paste0("Reading in data.set ", data.set, ": ", file, "\n"))
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
    
    
    ## Patients listed in the annotation file may have 0, 1, or >1 corresponding
    ## samples in the data file.  Multiple samples will be semicolon delimited.
    ## If there are multiple, take the mean across samples.
    ## (NB: since most of these are log expression values, this corresponds
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

    ## Combine all of the sample name maps into a single map
    map <- do.call("rbind", sample.name.mappings)
    rownames(map) <- NULL
    
    map <- unique(map)
    if(any(duplicated(map$patient.id))) {
        warning("Was not expecting any patient IDs to be duplicated\n")
    }
    
    ## The sample IDs are scattered across the multiple data type columns,
    ## put them in the single "sample.id" column
    anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)

    ## Subset the expression data sets to only have genes that are common
    ## to all data sets.
    all.genes <- lapply(data.sets, rownames)
    common.genes <- Reduce(intersect, all.genes)
    
    ## Subset all of the data sets to have only those genes shared
    ## across all data sets.
    for(data.set in names(data.sets)) {
        data.sets[[data.set]] <- data.sets[[data.set]][as.character(common.genes),]
    }
    
    ## Z-score each gene, separately for each data set.
    ## NB: scale normalizes the _columns_ of the matrix
    for(data.set in names(data.sets)) {
        data.sets[[data.set]] <- t(scale(t(data.sets[[data.set]]), center = TRUE, scale = TRUE))
    }
    
    ## Combine all of the expression matrices into a single matrix and
    ## ensure they are ordered consistently with the survival data.
    expr <- do.call("cbind", data.sets)
    common.samples <- intersect(colnames(expr), anno.tbl$sample.id)
    anno.tbl <- subset(anno.tbl, sample.id %in% common.samples)
    expr <- expr[, anno.tbl$sample.id]

    list(expr = expr, anno.tbl = anno.tbl)
}
    
## Read in the annotation file, which indicates the files
## holding the training data.
anno.tbl <- read.annotation.file()

lst <- process.expression.data.sets(anno.tbl)
anno.tbl <- lst[["anno.tbl"]]
expr <- lst[["expr"]]
expr <- expr[!unlist(apply(expr, 1, function(row) any(!is.finite(row)))),]

## Assemble the clinical data that are features in the model

clinical.features <- anno.tbl[,c("D_Gender", "D_ISS", "D_Age", "Study")]
rownames(clinical.features) <- anno.tbl$sample.id

## Let's impute gender using the y-chromosome genes
na.flag <- is.na(clinical.features$D_Gender)
if(any(na.flag)) {
    ## Look for strong associations between gender and expression of genes
    ## on the Y chromosome (these were stored during the training phase
    ## in y.genes)
    gender.pvals <- ldply(intersect(y.genes, rownames(expr)), .parallel = FALSE, .fun = function(gene) c(gene, wilcox.test(expr[gene,] ~ anno.tbl$D_Gender)$p.val))
    y.genes <- gender.pvals$V1[gender.pvals$V2 < 10^-10]

    tmp <- cbind(as.data.frame(t(expr[intersect(y.genes,rownames(expr)),,drop=F])), gender = anno.tbl$D_Gender)
    is.male <- tmp$gender == "Male"
    tmp$gender <- as.factor(tmp$gender)
    tmp$ng <- NA
    tmp$ng[is.male] <- 0
    tmp$ng[!is.male] <- 1
    ca <- preProcess(tmp[, !(colnames(tmp) == "gender")], method="medianImpute")
    predictions <- predict(ca, tmp)
    clinical.features$D_Gender_Imputed <- unlist(lapply(predictions$ng, function(x) ifelse(x == 0, "Male", "Female")))
    cat("Imputed Gender\n")
    table(clinical.features$D_Gender, clinical.features$D_Gender_Imputed)
    cat(paste0("Imputing D_Gender for ", length(which(na.flag)), " of ", nrow(clinical.features), " samples\n"))
    clinical.features$D_Gender[na.flag] <- clinical.features$D_Gender_Imputed[na.flag]
}
clinical.features$D_Gender <- factor(clinical.features$D_Gender)

## Impute ISS
## No, don't do imputation.  The study GSE15695 has a majority of patients
## that are ISS=3.  This means we will impute the NAs as 3, which will
## lead us to call those patients high risk.  Instead, "impute" them
## to have ISS = 1, which means we will just rely on EMC-92.
clinical.features$D_ISS <- as.numeric(clinical.features$D_ISS)
cat(paste0("ISS: ", paste(unique(clinical.features$D_ISS), collapse=", "), "\n"))
na.flag <- is.na(clinical.features$D_ISS)
if(any(na.flag)) {
    cat(paste0("Imputing D_ISS for ", length(which(na.flag)), " of ", nrow(clinical.features), " samples\n"))
    clinical.features$D_ISS[na.flag] <- 1
}
clinical.features$D_ISS <- factor(clinical.features$D_ISS)

## Impute age
na.flag <- is.na(clinical.features$D_Age)
if(any(na.flag)) {
    cat(paste0("Imputing D_Age for ", length(which(na.flag)), " of ", nrow(clinical.features), " samples\n"))
    ## Assign it the median age in its data set
    for(study in unique(clinical.features$Study)) {
        study.flag <- clinical.features$Study == study
        if(any(study.flag & na.flag)) {
            median.age <- median(as.numeric(clinical.features$D_Age[study.flag]), na.rm=TRUE)
            cat(paste0("Imputing D_Age for study ", study, " as median age = ", median.age, "\n"))
            clinical.features$D_Age[na.flag & study.flag] <- median.age
        }
    }
}

clinical.features <- clinical.features[, c("D_ISS", "D_Gender", "D_Age")]
texpr <- t(expr[intersect(target.gene.entrez,rownames(expr)),,drop=F])
common.samples <- intersect(rownames(texpr), rownames(clinical.features))
texpr <- texpr[common.samples,]
clinical.features <- clinical.features[common.samples,]
data <- cbind(texpr, clinical.features)

cat(paste0("any(is.na(data)) = ", any(is.na(data)), "\n"))

cat("Calling prediction\n")
predicted <- predict(trained.rf, newdata = data)$predicted

predictions.tbl <- data.frame(predictionscore=predicted, patient=rownames(data))
pfs.threshold.days <- (365.25/12)*18

## The rfsrc prediction output is a "mortality".
## I don't know how this relates to high-risk, per se.  I believe we would need to train a threshold.
## However, I believe higher mortality correlates with higher-risk (i.e., lower PFS)
predictions.tbl$highriskflag <- unlist(lapply(predictions.tbl$predictionscore, function(x) as.numeric(x) < pfs.threshold.days))
## predictions.tbl$predictionscore <- 1/as.numeric(predictions.tbl$predictionscore)
## Changed on 7/5/17
predictions.tbl$predictionscore <- as.numeric(predictions.tbl$predictionscore)

print(head(predictions.tbl))

dir.create(output.dir)

## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")
predictions.tbl <- predictions.tbl[ ,c("patient", "predictionscore", "highriskflag")]

## Merge in the study name, which is required in the prediction output.
predictions.tbl <- merge(predictions.tbl, anno.tbl[,c("Patient", "Study")], by.x = "patient", by.y = "Patient")
predictions.tbl <- predictions.tbl[,c("Study", "patient", "predictionscore", "highriskflag")]
colnames(predictions.tbl) <- c("study", "patient", "predictionscore", "highriskflag")

write.table(file=output.file, predictions.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")

