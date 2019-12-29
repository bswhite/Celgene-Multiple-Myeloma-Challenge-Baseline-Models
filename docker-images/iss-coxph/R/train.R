#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("maxstat"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix.utils"))

option_list <- list(
                    make_option(c("--sc"), action="store",
                                default=NULL,
                                help="String indicating which sub-challenge to train (sc1, sc2, or sc3)"),
                    make_option(c("--challenge-data-table-synId"), action="store",
                                default="syn9763946",
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

synapseLogin()

anno.file <- paste0(sc.string, "_Training_ClinAnnotations.csv")
docker.image.base <- "./"

data.dir <- paste0(docker.image.base, "../download/")
model.state.metadata.file <- paste0(docker.image.base, "/", sc.string, "_model-state-metadata.Rd")

dir.create(data.dir)

## Read in utility files.
common.src.dir <- paste0(docker.image.base, "/R/")
common.src.dir <- "../../publishedClassifiers/common/"
source(paste0(common.src.dir, "impute.R"))
source(paste0(common.src.dir, "utilities.R"))

## Read in the annotation data
anno.tbl <- read.annotation.data(anno.file, data.dir, download.data = download.data, synId = challenge.data.table.synId)

## Impute Age, Gender, and ISS as required.
anno.tbl <- impute.missing.clinical.data(anno.tbl)
colnames(anno.tbl)[colnames(anno.tbl) == "Patient"] <- "sample.id"

## Merge the model-based features with the clinical variables
clinical.covariates <- c("D_ISS")

y <- prepare.covariates(anno.tbl, clinical.covariates, include.events = TRUE)

res <- fit.and.threshold.coxph.model(y, clinical.covariates)
y <- res[["y"]]
model <- res[["fitted.model"]]
threshold <- res[["threshold"]]

if(FALSE) {
    y$high.risk <- unlist(apply(y[, c("event", "time")], 1,
                                function(row) ifelse((row[1] == 1) && (row[2] < 18*(365.25/12)), "high", "low")))
    y$inferred.high.risk <- y$cph > threshold
    fisher.test(y$high.risk, y$inferred.high.risk)

    library(ggplot2)
    ## library(ggbeeswarm)
    g <- ggplot(data = y, aes(x = high.risk, y = cph))
    g <- g + geom_boxplot()
    ## g <- g + geom_beeswarm()
}

## Save the trained state of the model, which consistts of
## (1) "threshold" -- the threshold for discretizing the continuous prediction into high vs low risk.
## (2) "model" -- the coxph model, which can be executed via 'survfit(model, newdata = y)'

model.state <- list("threshold" = threshold, "model" = model)

save(model.state, file=model.state.metadata.file)

cat("Successfully trained model\n")
