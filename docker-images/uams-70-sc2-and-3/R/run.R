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
Run a model for the sub-Challenge specified by the sc parameter.
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

## Read in the model implementation.
model.src.dir <- paste0(docker.image.base, "/R/")
if(download.data) {
    model.src.dir <- "../../publishedClassifiers/uams-70/"
}
source(paste0(common.src.dir, "classifier-base.R"))
source(paste0(model.src.dir, "uams-70.R"))

## Read in the model state
## It is assumed to be a list with two entries:
## "threshold"--the threshold for our classifier that was optimized on training data.
## "mapping"--the probe to entry mapping
## "model"--the (trained) coxph model
load(model.state.metadata.file)
threshold <- unname(model.state[["threshold"]])
model.entrezs <- model.state[["mapping"]]

genes.in.validation.data <- NULL

## Specifying download.data = TRUE and the synId of the sequestered data allows us to access
## that sequestered data and to avoid the need for Docker.
ret <- read.expression.and.annotation.data(anno.file, data.dir, download.data = download.data, synId = challenge.data.table.synId)
data.set.cols <- ret[["data.set.cols"]]
data.sets <- ret[["data.sets"]]
anno.tbl <- ret[["anno.tbl"]]

ret <- process.expression.and.annotation.data(anno.tbl, data.sets, data.set.cols, mapping = model.entrezs,
                                              genes.in.validation.data = genes.in.validation.data)
expr <- ret[["expr"]]
anno.tbl <- ret[["anno.tbl"]]

## Execute the UAMS-70 model using the trained features.
X <- as.matrix(expr)
y <- NULL
eset <- ExpressionSet(as.matrix(X))
coefficients <- uams.70.get.coefficients()
res <- uams.70.gene(eset, mapping = model.entrezs, coefficients = coefficients, anno = y, threshold = threshold)$res

## Assemble the predictions
predictions.tbl <- res[,c("ID", "raw.score", "high.risk")]
colnames(predictions.tbl) <- c("patient", "predictionscore", "highriskflag")

## Merge in the study name, which is required in the prediction output.
predictions.tbl <- merge(predictions.tbl, anno.tbl[,c("Patient", "Study")], by.x = "patient", by.y = "Patient")
predictions.tbl <- predictions.tbl[,c("Study", "patient", "predictionscore", "highriskflag")]
colnames(predictions.tbl) <- c("study", "patient", "predictionscore", "highriskflag")

## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")
dir.create(output.dir)
write.table(file=output.file, predictions.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")


