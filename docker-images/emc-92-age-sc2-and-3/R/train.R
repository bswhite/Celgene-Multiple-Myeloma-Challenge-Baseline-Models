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

## Read in the model implementation.
model.src.dir <- paste0(docker.image.base, "/R/")
model.src.dir <- "../../publishedClassifiers/emc-92/"
source(paste0(common.src.dir, "classifier-base.R"))
source(paste0(model.src.dir, "emc-92.R"))

## Get the Entrez genes in the model.  This is defined in the model implementation.
model.entrezs <- emc.92.probe.to.entrez.mapping()
model.entrezs <- na.omit(model.entrezs)

## Read in the list of genes that are present in the validation data.
## This will ensure that we train our model on genes that are guaranteed to be present in
## all data sets (training and validation).
genes.in.validation.data <- get.genes.in.validation.data(data.dir)

ret <- read.expression.and.annotation.data(anno.file, data.dir, download.data = download.data, synId = challenge.data.table.synId)
data.set.cols <- ret[["data.set.cols"]]
data.sets <- ret[["data.sets"]]
anno.tbl <- ret[["anno.tbl"]]

ret <- process.expression.and.annotation.data(anno.tbl, data.sets, data.set.cols, mapping = model.entrezs,
                                              genes.in.validation.data = genes.in.validation.data)
expr <- ret[["expr"]]
anno.tbl <- ret[["anno.tbl"]]

## Execute the EMC-92 model to define model-based features.
## Use a dummy threshold: we just want the continuous scores, not the thresholded/dichotomized calls
X <- as.matrix(expr)
y <- NULL
eset <- ExpressionSet(as.matrix(X))
coefficients <- emc.92.get.coefficients()
res <- emc.92.gene(eset, mapping = model.entrezs, coefficients = coefficients, anno = y, threshold = 0)$res

## Merge the model-based features with the clinical variables
clinical.covariates <- c("D_Age")
y <- anno.tbl[, c("sample.id", "D_PFS", "D_PFS_FLAG", clinical.covariates)]
## NB: Hovon has all Age = NA.  We will exclude that study in training.
y <- na.omit(y)
colnames(y) <- c("ID", "time", "event", clinical.covariates)
y <- merge(y, res[, c("ID", "raw.score")], by = "ID")
if("D_ISS" %in% colnames(y)) {
    y$D_ISS <- unlist(lapply(y$D_ISS,
                             function(x) {
                                 if(!(x %in% c(1,2,3))) { stop(paste0("Was not expecting ISS = ", x, "\n")) }
                                 c("I", "II", "III")[x]
                             }))
    y$D_ISS <- factor(y$D_ISS)
                      
}

## Fit a Cox proportional hazards model to the model-based feature and the clinical covariates
surv.formula <- as.formula(paste("Surv(time, event)", "~",
                                 paste(c("raw.score", clinical.covariates), collapse = " + "), sep = " "))
## surv.formula <- as.formula(paste("Surv(time, event)", "~", "raw.score", sep = " "))
cph.fit <- coxph(surv.formula, data = y)

pred <- unname(predict(cph.fit, newdata = y, type = "risk"))
y$cph <- pred

## Calculate a threshold for the Cox proportional hazards model: using the t = 0 
suppressPackageStartupMessages(library(maxstat))
mt.formula <- as.formula(paste("Surv(time, event)", "~", "cph", sep = " "))
mt <- maxstat.test(mt.formula, data=data.frame(y), smethod="LogRank", pmethod="none")
threshold <- as.numeric(mt$estimate)
model <- cph.fit

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
## (2) "mapping" -- a mapping from probe ids to entrez ids.  We will need this mapping during the
## validation phase where we don't have access to the internet.
## (3) "model" -- the coxph model, which can be executed via 'survfit(model, newdata = y)'

model.state <- list("threshold" = threshold, "mapping" = model.entrezs, "model" = model)

save(model.state, file=model.state.metadata.file)

cat("Successfully trained model\n")
