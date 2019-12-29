#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("maxstat"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("randomForest"))

synapseLogin()

## docker build -t docker.synapse.org/syn7237116/expr-prolif-cyto:version1 .
## docker images
## docker login docker.synapse.org
## docker push docker.synapse.org/syn7237116/expr-prolif-cyto:version1

base.dir <- "./"

download.dir <- paste0(base.dir, "download/")
dir.create(download.dir)

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

translate.ensg.to.entrez <- function(ids) {
  # Get a mapping from Ensembl id to entrez id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
              filters = 'ensembl_gene_id', 
              values = ids, 
              mart = ensembl)
  names(bm) <- c("ID", "TRANSLATION")
  bm <- bm[!(bm$TRANSLATION %in% c("")),]
  bm
}

translate.sym.to.entrez <- function(ids) {
  # Get a mapping from Ensembl id to entrez id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), 
              filters = 'hgnc_symbol', 
              values = ids, 
              mart = ensembl)
  names(bm) <- c("ID", "TRANSLATION")
  bm <- bm[!(bm$TRANSLATION %in% c("")),]
  bm
}

read.annotation.file <- function() {
    chal_data_table <- synTableQuery('select id,name from syn9763946')
    chal_data_df <- chal_data_table@values
    chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc2_Training_ClinAnnotations.csv"))
    sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=download.dir)})
    training.annotation.file <- paste0(download.dir, "sc2_Training_ClinAnnotations.csv")
    anno.tbl <- read.table(training.annotation.file, sep=",", header=TRUE)
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
    chal_data_table <- synTableQuery('select id,name from syn9763946')
    chal_data_df <- chal_data_table@values
    chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
    ## Download the data and put it in the base directory of the docker image,
    ## which is where it would be mounted in the challenge.
    sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=download.dir)})
    
    ## Read in each of the data sets and create:
    ## 1) A list holding each of the data sets.
    ## 2) A list with an entry for each data set holding a mapping from the
    ##    sample name in the clinical annotation column and the sample name
    ##    in the expression data set.
    data.sets <- list()
    sample.name.mappings <- list()
    
    for(col.indx in 1:length(data.set.cols)) {
        for(data.set in na.omit(unique(anno.tbl[,data.set.cols[col.indx]]))) {
            file <- paste0(download.dir, data.set)
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

    ## Translate Ensembl genes in the MMRF data set to Entrez.
    mmrf.data.set <- names(data.sets)[grepl(names(data.sets), pattern="MMRF")]
    
    ensg.to.entrez.map <- translate.ensg.to.entrez(rownames(data.sets[[mmrf.data.set]]))
    ensg.to.entrez.map <- na.omit(ensg.to.entrez.map)

    ## Subset the expression data sets to only have the genes in the UAMS-5
    ## data set (and that are in common between all data sets).
    
    all.genes <- lapply(data.sets, rownames)
    all.genes[[mmrf.data.set]] <- ensg.to.entrez.map$TRANSLATION[ensg.to.entrez.map$ID %in% all.genes[[mmrf.data.set]]]
    
    common.genes <- Reduce(intersect, all.genes)
    
    ## Finally, translate the MMRF Ensembl genes to Entrez
    ensg.to.entrez.map <- subset(ensg.to.entrez.map, TRANSLATION %in% common.genes)
    if(any(duplicated(ensg.to.entrez.map$ID))) {
        warning("Ambiguous mapping between Ensembl ID an Entrez TRANSLATION\n")
        ensg.to.entrez.map <- ensg.to.entrez.map[!duplicated(ensg.to.entrez.map$ID, fromLast=TRUE) & !duplicated(ensg.to.entrez.map$ID, fromLast=FALSE),]
    }
    common.genes <- unique(ensg.to.entrez.map$TRANSLATION)
    
    data.sets[[mmrf.data.set]] <- data.sets[[mmrf.data.set]][rownames(data.sets[[mmrf.data.set]]) %in% ensg.to.entrez.map$ID,]
    data.sets[[mmrf.data.set]] <- data.sets[[mmrf.data.set]][ensg.to.entrez.map$ID,]
    rownames(data.sets[[mmrf.data.set]]) <- ensg.to.entrez.map$TRANSLATION
    
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
    ## ensure they are ordered consistently with the survival data.
    expr <- do.call("cbind", data.sets)
    expr <- expr[, anno.tbl$sample.id]

    list(expr = expr, anno.tbl = anno.tbl, ensg.to.entrez.map = ensg.to.entrez.map)
}
    
## Read in the annotation file, which indicates the files
## holding the training data.
anno.tbl <- read.annotation.file()

lst <- process.expression.data.sets(anno.tbl)
anno.tbl <- lst[["anno.tbl"]]
expr <- lst[["expr"]]
expr <- expr[!unlist(apply(expr, 1, function(row) any(!is.finite(row)))),]

## Train a RF model using age, gender, iss, and genes that are targets of known myeloma cytogenetic
## markers (IGH and MYC translocations and TP53 deletions).  We will approximate these cytogenetic
## abnormalities by looking at the expression of the IGH/MYC targets and of TP53.

## 1. Impute any missing Gender
## Get genes on the Y chromosomes
library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes = c("chromosome_name", "entrezgene", "hgnc_symbol"), filters = "chromosome_name", values = c("Y"), mart = mart)
y.genes <- results$entrezgene

## Look for strong associations between expression and gender
gender.pvals <- ldply(intersect(y.genes, rownames(expr)), .parallel = FALSE, .fun = function(gene) c(gene, wilcox.test(expr[gene,] ~ anno.tbl$D_Gender)$p.val))
y.genes <- gender.pvals$V1[gender.pvals$V2 < 10^-10]

tmp <- cbind(as.data.frame(t(expr[intersect(y.genes,rownames(expr)),,drop=F])), gender = anno.tbl$D_Gender)
tmp$gender <- as.factor(tmp$gender)

## suppressPackageStartupMessages(library("missForest"))
## mf <- missForest(tmp, replace=TRUE, mtry = ncol(tmp)-1)
## anno.tbl$D_Gender_Imputed <- mf$ximp[anno.tbl$sample.id,"gender"]
## table(anno.tbl$D_Gender, anno.tbl$D_Gender_Imputed)

## fit <- glm(gender ~ . , family = binomial("logit"), data = tmp)
## probabilities <- predict(fit, tmp[, !(colnames(tmp) == "gender")], type="response")
## predictions <- probabilities
## predictions[probabilities < 0.5] <- levels(tmp$gender)[1]
## predictions[probabilities > 0.5] <- levels(tmp$gender)[2]
## predictions <- factor(predictions, levels = levels(tmp$gender))
## table(predictions, tmp$gender)

## suppressPackageStartupMessages(library("caret"))
## ca <- preProcess(tmp, method="knnImpute")
## predictions <- predict(ca, tmp)
## anno.tbl$D_Gender_Imputed <- predictions$gender
## table(anno.tbl$D_Gender, anno.tbl$D_Gender_Imputed)
## anno.tbl$D_Gender[is.na(anno.tbl$D_Gender)] <- predictions[is.na(anno.tbl$D_Gender),"gender"]

## NSD2 = MMSET/WHSC1
target.gene.symbols <- c("MYC", "TP53", "FGFR3", "WHSC1", "MMSET", "CCND3", "CCND1", "MAF", "MAFB", "NSD2")

proliferation.gene.symbols <- c("KIF4A", "TPX2", "DLGAP5", "NUF2", "SPAG5", "TOP2A", "ASPM", "BUB1", "POLQ", "SGO2", "ANLN", "KIF23", "KIF21B", "CENPF", "ESPL1", "TACC3", "MKI67", "DTL", "PLK4", "NCAPH", "KIF14", "BRCA1", "MCM10", "CENPE", "CEP55", "MCM4", "NEK2", "PHF19", "CDCA2", "KIF15", "KIF11", "CDT1", "FANCI", "TK1", "KNL1", "EXO1", "E2F2", "CCDC15", "BIRC5", "TTK", "EZH2", "HJURP", "PRR11", "NCAPG", "ZWINT", "AURKA", "FOXM1", "NUSAP1", "UBE2T", "SMC2", "HMMR", "GPR63", "KIF20A", "TRIP13", "DEPDC1", "CDCA8", "KIF2C", "MCM2", "GAS2L3", "PROSER3", "E2F8", "NCAPG2", "STIL", "SKA1", "CDC45", "UBE2C", "CDC20", "DIAPH3", "WHSC1", "NDE1", "RAD18", "PKP2", "CENPU", "UBAP2L", "POLD3", "POLD1", "STEAP1", "ATAD2", "MTFR2", "DSCC1", "HMGN5", "SPDL1", "POLA1", "SMC4", "KIF18B", "CCNB2", "GINS4", "DEPDC1B", "CDKN3", "NDC80", "MELK", "SPC24", "CENPN", "TROAP", "HELLS", "PCYOX1L", "CDK1", "SGO1", "DARS2", "CENPW")

ribosomal.gene.symbols <- c("REEP5", "RPL18", "PLPP5", "SAT2", "KIAA1191", "RPL3", "RPS2", "TPT1", "DNAJB9", "RPL24", "ITM2B", "EEF1A1", "MTCH1", "CSNK1G3", "LETMD1", "TSEN34", "TM9SF1", "FAU", "TM9SF2", "RPS17", "RPL12", "EIF2A", "TMEM243", "RPS19", "RPL30", "ATP5L", "NSA2", "TMEM258", "GPR160", "VPS51", "TNFRSF14", "PAM", "SURF4", "NANS", "RPL10A", "CREBL2", "PHF1", "RPS16", "RPL35A", "RPL9", "RPL41", "STAP1", "GABARAP", "RPL34", "NDFIP2", "RPS25", "TSPAN31", "HSCB", "RPL31", "RPS6", "TXNDC15", "PHB2", "RPS4X", "EIF3L", "UBE2D2", "ISCU", "RPL8", "PPP1CC", "RPL14", "SLC25A6", "CAMLG", "C2orf88", "RPL27", "SPCS2", "TLR10", "OSBPL10", "RPL22", "BTD", "CDV3", "TBCA", "SEL1L3", "B4GAT1", "CRYL1", "AP3S1", "RPS3", "TMBIM4", "KDELR1", "SLC7A7", "SRPRB", "ZFAND1", "LYRM1", "RPL15", "TAPBPL", "SERP1", "CCNG1", "RPL7", "AMIGO2", "MRFAP1", "RPLP1", "SLC31A2", "SAR1B", "COX7C", "RPL36AL", "ZDHHC24", "TESK2", "PRCP", "RNF11", "SEC24A", "PTTG1IP")

tmp <- translate.sym.to.entrez(target.gene.symbols)
target.gene.entrez <- tmp$TRANSLATION

tmp <- translate.sym.to.entrez(proliferation.gene.symbols)
proliferation.gene.entrez <- tmp$TRANSLATION

tmp <- translate.sym.to.entrez(ribosomal.gene.symbols)
ribosomal.gene.entrez <- tmp$TRANSLATION

all.entrezs <- c(target.gene.entrez, proliferation.gene.entrez, ribosomal.gene.entrez)

## 18 months
pfs.threshold.days <- (365.25/12)*18
anno.tbl$pfsBasedHighRisk <- unlist(apply(anno.tbl[, c("D_PFS", "D_PFS_FLAG")], 1, function(row) ifelse(any(is.na(row)), NA, ifelse((row[1] < pfs.threshold.days) & (row[2] == 1), TRUE, FALSE))))
clinical.features <- anno.tbl[,c("D_Gender", "D_ISS", "D_Age", "pfsBasedHighRisk")]
rownames(clinical.features) <- anno.tbl$sample.id
## texpr <- t(expr[intersect(target.gene.entrez,rownames(expr)),,drop=F])
ex <- expr[intersect(target.gene.entrez, rownames(expr)),,drop=F]
prolif.avg <- colMeans(expr[intersect(proliferation.gene.entrez, rownames(expr)),,drop=F])
ribo.avg <- colMeans(expr[intersect(ribosomal.gene.entrez, rownames(expr)),,drop=F])
ex <- rbind(ex, prolif.avg, ribo.avg)
texpr <- t(ex)

common.samples <- intersect(rownames(texpr), rownames(clinical.features))
texpr <- texpr[common.samples,]
clinical.features <- clinical.features[common.samples,]
data <- cbind(texpr, clinical.features)
data <- na.omit(data)
data$D_ISS <- factor(data$D_ISS)
data$D_Gender <- factor(data$D_Gender)
data$pfsBasedHighRisk <- factor(data$pfsBasedHighRisk)

lbl <- data$pfsBasedHighRisk
min.factor <- min(table(lbl))
max.factor <- max(table(lbl))
num.factors <- length(unique(lbl))
sampsize <- rep(floor(0.7*min.factor), num.factors)

trained.rf <- randomForest(x=data[, !(colnames(data) == "pfsBasedHighRisk")], y=lbl, keep.forest=TRUE, proximity=FALSE, strata=lbl, sampsize=sampsize, importance=TRUE, replace=FALSE)

## trained.rf <- train(pfsBasedHighRisk ~ ., data = data, method='rf', strata=lbl, sampsize=sampsize, replace=FALSE)
cat(paste0("Class of trained.rf = ", class(trained.rf), "\n"))

model.state <- list("model" = trained.rf, "target.gene.entrez" = target.gene.entrez, "proliferation.gene.entrez" = proliferation.gene.entrez, "ribosomal.gene.entrez" = ribosomal.gene.entrez, "y.genes" = y.genes)

model.state.metadata.file <- paste0(base.dir, "/", "model-state-metadata.Rd")
save(model.state, file=model.state.metadata.file)

cat("Successfully trained model\n")

