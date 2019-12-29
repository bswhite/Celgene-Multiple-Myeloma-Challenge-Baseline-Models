#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix.utils"))

option_list <- list(
                    make_option(c("--sc"), action="store",
                                default=NULL,
                                help="String indicating which sub-challenge to train (sc1, sc2, or sc3)"),
                    make_option(c("--challenge-data-table-synId"), action="store",
                                default=NULL,
                                help="The Synapse ID of the table enumerating the challenge data (or NULL if data should instead by accessed locally).  Default: %default%")
    )

descr <- "\
Train a model for the sub-Challenge specified by the sc parameter.
"

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if ( length(arguments$args) != 0 ) {
  print_help(parser)
  q(status=1)
}

## Collect input parameters
sc.string <- opt$sc

if(!(sc.string %in% c("sc1", "sc2", "sc3"))) {
    stop(paste0("Unexpected arg: ", sc.string, "\n"))
}

challenge.data.table.synId <- opt$`challenge-data-table-synId`
download.data <- !is.null(challenge.data.table.synId)

if(download.data) {
    suppressPackageStartupMessages(library("synapseClient"))
    suppressPackageStartupMessages(library("maxstat"))

    synapseLogin()
}

anno.file <- paste0(sc.string, "_Validation_ClinAnnotations.csv")

docker.image.base <- "./"

data.dir <- paste0(docker.image.base, "test-data/")
output.dir <- paste0(docker.image.base, "output/")
model.state.metadata.file <- paste0(docker.image.base, "/", sc.string, "_model-state-metadata.Rd")

if(download.data) {
    anno.file <- paste0(sc.string, "_INTERNALVALIDATION_clinAnnotation.csv")
    dir.create(data.dir)
}

## Read in utility files.
common.src.dir <- paste0(docker.image.base, "/R/")
if(download.data) {
    common.src.dir <- "../../publishedClassifiers/common/"
}
source(paste0(common.src.dir, "impute.R"))
source(paste0(common.src.dir, "utilities.R"))

## Read in the model state
## It is assumed to be a list with two entries:
## "threshold"--the threshold for our classifier that was optimized on training data.
## "model"--the (trained) coxph model
load(model.state.metadata.file)
threshold <- unname(model.state[["threshold"]])
model <- model.state[["model"]]
print(summary(model))

## Read in the annotation data
anno.tbl <- read.annotation.data(anno.file, data.dir, download.data = download.data, synId = challenge.data.table.synId)
colnames(anno.tbl)[colnames(anno.tbl) == "Patient"] <- "sample.id"

## Impute Age, Gender, and ISS as required.
anno.tbl <- impute.missing.clinical.data(anno.tbl)

## Merge the model-based features with the clinical variables
clinical.covariates <- c("D_ISS", "D_Age")

y <- prepare.covariates(anno.tbl, clinical.covariates)

pred <- unname(predict(model, newdata = y, type = "risk"))
y$raw.score <- pred
y$high.risk <- y$raw.score > threshold

predictions.tbl <- y[,c("ID", "raw.score", "high.risk")]
colnames(predictions.tbl) <- c("patient", "predictionscore", "highriskflag")

## Merge in the study name, which is required in the prediction output.
predictions.tbl <- merge(predictions.tbl, anno.tbl[,c("sample.id", "Study")], by.x = "patient", by.y = "sample.id")
predictions.tbl <- predictions.tbl[,c("Study", "patient", "predictionscore", "highriskflag")]
colnames(predictions.tbl) <- c("study", "patient", "predictionscore", "highriskflag")

dir.create(output.dir)

## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")

write.table(file=output.file, predictions.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")


