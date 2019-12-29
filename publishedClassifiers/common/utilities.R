## Utility functions for assessing challenge data

get.challenge.file.names <- function(synId = "syn9763946") {
    chal_data_table <- synTableQuery(paste0("select id,name from ", synId))
    chal_data_df <- chal_data_table@values
    chal_data_df
}

download.challenge.data <- function(chal_data_df, download.dir) {
    sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=download.dir)})
}

read.challenge.data.sets <- function(anno.tbl, data.set.cols, download.dir) {

    ## Get each (unique) data set file listed in any of the data.set.cols columns in the
    ## annotation file.
    names(data.set.cols) <- data.set.cols
    df <- ldply(data.set.cols, .fun = function(col) data.frame(na.omit(unique(as.character(anno.tbl[, col])))))
    colnames(df) <- c("data.set.col.id", "data.set.file")
    df <- unique(df)
    
    ## Read in each file.
    data.set.files <- as.character(unique(df$data.set.file))
    names(data.set.files) <- data.set.files

    data.sets <- llply(data.set.files,
                       .fun = function(data.set) {
                           cat(paste0("Reading in data.set ", data.set, "\n"))
                           file <- paste0(download.dir, data.set)
                           tbl <- as.data.frame(fread(file, header=TRUE))

                           ## The first column is the gene identifier--make it the row
                           ## name and drop it from the table
                           rownames(tbl) <- tbl[,1]
                           tbl <- tbl[,-1]
                           tbl
                       })
    data.sets
}

create.sample.name.mapping <- function(anno.tbl, data.set.cols) {

    ## Get each (unique) data set file listed in any of the data.set.cols columns in the
    ## annotation file.
    names(data.set.cols) <- data.set.cols
    df <- ldply(data.set.cols, .fun = function(col) data.frame(na.omit(unique(as.character(anno.tbl[, col])))))
    colnames(df) <- c("data.set.col.id", "data.set.file")
    df <- unique(df)

    if(any(duplicated(df$data.set.file))) {
        stop("Was not expecting any files to be duplicated across data set cols\n")
    }

    indices <- 1:nrow(df)
    names(indices) <- df$data.set.file

    sample.name.mappings <- llply(indices,
                                  .fun = function(indx) {
                                      ## Extract the sample name mappings from the annotation file
                                      col <- as.character(df$data.set.col.id[indx])
                                      data.set <- as.character(df$data.set.file[indx])
                                      
                                      flag <- !is.na(anno.tbl[,col]) & (as.character(anno.tbl[,col]) == data.set)

                                      ## The mapping between the Patient/row of the clinical annotation file and
                                      ## the identifier in each data file is provided by the corresponding ID
                                      ## columns of the clinical annotation file.
                                      ## These have the same names as the data set columns, but with "SamplId"
                                      ## appended.
                                      patient.id.col <- paste0(col, "SamplId")
                                      map <- anno.tbl[flag,c("Patient", patient.id.col)]
                                      colnames(map) <- c("patient.id", "sample.id")
                                      map
                                  })
    sample.name.mappings
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

## A Patient may be associated with multiple SamplIds, as occurs when multiple samples from a
## patient were assayed. The multiple SamplIds will be semicolon-delimited. For example,
## Patient "MMRF_1024" has three SamplIds listed in the RNASeq_geneLevelExpFileSamplId
## column: "MMRF_1024_1_BM;MMRF_1024_2_BM;MMRF_1024_3_BM".
## These correspond to three different bone marrow specimens from the Patient, all of which
## have data in the expression file.

## Patients listed in the annotation file may have 0, 1, or >1 corresponding
## samples in the data file.  If there are multiple, take the mean across
## samples (NB: since most of these are log expression values, this corresponds
## to the geometric mean in real space).
## Also, subset both the mapping tables and the data to have the
## intersection of patients in both.
handle.multiple.samples.per.patient <- function(data.sets, sample.name.mappings) {

    ret.data.sets <- list()
    ret.sample.name.mappings <- list()
    
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
        
        ret.data.sets[[data.set]] <- tbl
        ret.sample.name.mappings[[data.set]] <- map
    }
    list(data.sets = ret.data.sets, sample.name.mappings = ret.sample.name.mappings)
}

## Merge the sample-to-patient mapping with the patient-level annotation data.
merge.sample.mappings.with.annotations <- function(anno.tbl, sample.name.mappings) {
    map <- do.call("rbind", sample.name.mappings)
    rownames(map) <- NULL
    map <- unique(map)
    if(any(duplicated(map$patient.id))) {
        warning("Was not expecting any patient IDs to be duplicated\n")
    }
    
    ## The sample IDs are scattered across multiple columns, put them
    ## in the single "sample.id" column
    anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)
    anno.tbl
}

## Get the genes in the validation data
get.genes.in.validation.data <- function(download.dir) {

    synId <- "syn10792130"
    obj <- synGet(synId, downloadLocation=download.dir)
    genes <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE)
    genes$entrezID

}

read.annotation.data <- function(anno.file, data.dir, download.data = FALSE, synId = "syn9763946") {
    chal_data_df <- NULL

    ## Download the names of the Challenge data files.
    if(download.data) {
        chal_data_df <- get.challenge.file.names(synId = synId)

        chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern=anno.file))
        download.challenge.data(chal_data_df_anno, data.dir)
    }
    
    ## Read in the annotation file, which indicates the files
    ## holding the training data.
    anno.file <- paste0(data.dir, anno.file)
    anno.tbl <- read.table(anno.file, sep=",", header=TRUE, as.is=TRUE)
    anno.tbl
}

read.expression.and.annotation.data <- function(anno.file, data.dir, download.data = FALSE, synId = "syn9763946") {
    chal_data_df <- NULL

    ## Download the names of the Challenge data files.
    if(download.data) {
        chal_data_df <- get.challenge.file.names(synId = synId)

        chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern=anno.file))
        download.challenge.data(chal_data_df_anno, data.dir)
    }
    
    ## Read in the annotation file, which indicates the files
    ## holding the training data.
    anno.file <- paste0(data.dir, anno.file)
    anno.tbl <- read.table(anno.file, sep=",", header=TRUE, as.is=TRUE)

    ## Read in the clinical annotation file that lists the samples and what data
    ## are available for each in the columns:
    data.set.cols <- colnames(anno.tbl)[grepl(pattern="File", colnames(anno.tbl)) & !grepl(pattern="Sampl", colnames(anno.tbl))]

    ## Restrict to gene-level expression data sets
    expression.data.set.flag <- grepl(pattern="geneLevelExp", data.set.cols)
    data.set.cols <- data.set.cols[expression.data.set.flag]

    if(download.data) {
        ## Only download the training data that is pertinent to this subchallenge.
        all.data.sets <- unique(na.omit(as.vector(as.matrix(anno.tbl[,data.set.cols]))))
        chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
    
        ## Download the data and put it in the base directory of the docker image,
        ## which is where it would be mounted in the challenge.
        download.challenge.data(chal_data_df, data.dir)
    }
    
    ## Read in each of the data sets and create a list 'data.sets' holding each of them.
    data.sets <- read.challenge.data.sets(anno.tbl, data.set.cols, data.dir)

    list(data.set.cols = data.set.cols, data.sets = data.sets, anno.tbl = anno.tbl)
}

process.expression.and.annotation.data <- function(anno.tbl, data.sets, data.set.cols, mapping,
                                                   genes.in.validation.data = NULL) {
    ## Create a list 'sample.name.mappings' with an entry for each data set file holding a mapping from the
    ## sample name in the clinical annotation column and the sample name
    ## in the expression data set.
    sample.name.mappings <- create.sample.name.mapping(anno.tbl, data.set.cols)
    
    ## Multiple sample ids may map to the same patient.  If so, take the mean of the sample expression
    ## values as the value for the patient.  Update the data.sets and sample.name.mappings accordingly
    ret <- handle.multiple.samples.per.patient(data.sets, sample.name.mappings)
    data.sets <- ret[["data.sets"]]
    sample.name.mappings <- ret[["sample.name.mappings"]]

    ## Merge the sample-to-patient mapping with the patient-level annotation data.
    anno.tbl <- merge.sample.mappings.with.annotations(anno.tbl, sample.name.mappings)
    
    ## Subset the expression data sets to only have the genes in the UAMS-70
    ## data set (and that are in common between all data sets).
    
    all.genes <- lapply(data.sets, rownames)
    
    common.genes <- Reduce(intersect, all.genes)
    common.genes <- intersect(common.genes, mapping$INDEX)
    
    ## Also, intersect with genes present in the validation data.
    if(!is.null(genes.in.validation.data)) {
        common.genes <- intersect(common.genes, genes.in.validation.data)
    }
    
    ## Subset all of the data sets
    for(data.set in names(data.sets)) {
        data.sets[[data.set]] <- data.sets[[data.set]][as.character(common.genes),]
    }

    ## Z-score each gene, separately for each data set.
    ## NB: scale normalizes the _columns_ of the matrix
    for(data.set in names(data.sets)) {
        data.sets[[data.set]] <- t(scale(t(data.sets[[data.set]]), center = TRUE, scale = TRUE))
    }

    ## Combine all of the expression matrices into a single matrix and
    ## ensure they are ordered consistently with the annotation data.
    expr <- do.call("cbind", data.sets)
    common.samples <- unique(intersect(colnames(expr), anno.tbl$sample.id))
    anno.tbl <- anno.tbl[anno.tbl$sample.id %in% common.samples,]
    expr <- expr[, anno.tbl$sample.id]
    
    ## Impute ISS, age, and gender
    anno.tbl <- impute.missing.clinical.data(anno.tbl)

    list(anno.tbl = anno.tbl, expr = expr)
}



prepare.covariates <- function(anno.tbl, clinical.covariates, include.events = FALSE) {

    y <- NULL
    if(include.events) {
        y <- anno.tbl[, c("sample.id", "D_PFS", "D_PFS_FLAG", clinical.covariates)]
    } else {
        y <- anno.tbl[, c("sample.id", clinical.covariates)]
    }

    ## NB: Hovon has all Age = NA.  We will exclude that study in training.                                                                                                            y <- na.omit(y)

    if(include.events) {
        colnames(y) <- c("ID", "time", "event", clinical.covariates)
    } else {
        colnames(y) <- c("ID", clinical.covariates)
    }
    if("D_ISS" %in% colnames(y)) {
        y$D_ISS <- unlist(lapply(y$D_ISS,
                                 function(x) {
                                     if(!(x %in% c(1,2,3))) { stop(paste0("Was not expecting ISS = ", x, "\n")) }
                                     c("I", "II", "III")[x]
                                 }))
        y$D_ISS <- factor(y$D_ISS)
    }
    y
}

fit.and.threshold.coxph.model <- function(y, clinical.covariates) {

    ## Fit a Cox proportional hazards model
    surv.formula <- as.formula(paste("Surv(time, event)", "~",
                                     paste(clinical.covariates, collapse = " + "), sep = " "))
    cph.fit <- coxph(surv.formula, data = y)

    pred <- unname(predict(cph.fit, newdata = y, type = "risk"))
    y$cph <- pred

    ## Calculate a threshold for the Cox proportional hazards model                                                                                                                     
    suppressPackageStartupMessages(library(maxstat))
    mt.formula <- as.formula(paste("Surv(time, event)", "~", "cph", sep = " "))
    mt <- maxstat.test(mt.formula, data=data.frame(y), smethod="LogRank", pmethod="none")

    list(y = y, fitted.model = cph.fit, threshold = as.numeric(mt$estimate))

}


