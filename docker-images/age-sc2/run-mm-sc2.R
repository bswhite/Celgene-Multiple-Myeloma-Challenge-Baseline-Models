#!/usr/bin/env Rscript

## R library dependencies. These must be installed in the
## Docker image (i.e., specified in Dockerfile)
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))

## Read in the trained model state, which is assumed to be
## in the model.state variable.
model.state.metadata.file <- "/model-state-metadata.Rd"
load(model.state.metadata.file)

## model.state is assumed to be a list with a single entry:
## "age"--the age cutoff above which patients should be
##        considered high risk
age.cutoff <- as.numeric(model.state[["age"]])
cat(paste0("Calling any patient with age > ", age.cutoff, " high risk\n"))
cat(paste0("Using age as prediction score\n"))

## Read in the annotation file for Sub-Challenge 2
anno.filename <- "/test-data/sc2_Validation_ClinAnnotations.csv"
cat(paste0("\nReading in annotation file ", anno.filename, "\n"))
anno.tbl <- fread(anno.filename)

## Use age is our prediction score--i.e., assume that higher age 
## is associated with higher risk of MM.
predictions <- as.data.frame(anno.tbl[, c("Study", "Patient", "D_Age")])

## Some patient ages may not be specified.
## "Impute" by setting to mean age of all data.
predictions[is.na(predictions$D_Age), "D_Age"] <-
    mean(as.numeric(predictions$D_Age), na.rm=TRUE)

## Binarize our continuous score (i.e., age) such that any above
## the "trained" threshold are called high risk
predictions$highriskflag <- FALSE

## This is a gratutious use of the plyr library
exceed.cutoff <- unlist(llply(as.numeric(predictions$D_Age), .fun = function(x) x > age.cutoff))
predictions$highriskflag[exceed.cutoff] <- TRUE

## Prediction table columns need to be named study, patient,
## predictionscore (continuous), and highriskflag (binarized)
old.col.names <- c("Study", "Patient", "D_Age", "highriskflag")
new.col.names <- c("study", "patient", "predictionscore", "highriskflag")
predictions <- predictions[, old.col.names]
colnames(predictions) <- new.col.names

## Output the results into /output/predictions.tsv
output.file <- "/output/predictions.tsv"
write.table(file=output.file, predictions, row.names=F, quote=F, sep="\t")

cat("Successfully wrote predictions.\n")

