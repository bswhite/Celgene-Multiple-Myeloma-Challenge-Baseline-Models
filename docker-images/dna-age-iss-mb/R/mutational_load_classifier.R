library(vcfR)
library(scales)

args = commandArgs(trailingOnly=TRUE)

testdir = '/test-data'
outputdir = '/output'

test.classifier <- function(input.file=file.path(testdir,args[1]),
mutationburden=64.51,
mutationFileField=c('WES_mutationFileMutect','WES_mutationFileStrelkaSNV','RNASeq_mutationFileMutect'),output.file=file.path(outputdir,'predictions.tsv')) {
  clinical_file <- read.csv(file=input.file,header=T,as.is=T,check.names=F)
  clinical_file$D_ISS[is.na(clinical_file$D_ISS)] <- median(as.numeric(clinical_file$D_ISS),na.rm=T);
  clinical_file$D_Age <- as.numeric(clinical_file$D_Age)
  clinical_file$D_Age[is.na(clinical_file$D_Age)] <- median(as.numeric(clinical_file$D_Age),na.rm=T);

  getMutationCount <- function(X,testdir) {
    tryCatch({
      X <- gsub("FILTERED.","",X)
      curr <- read.vcfR(file.path(testdir,X),verbose=FALSE);
      100*length(which(grepl(pattern='missense',curr@fix[,8]) & getFIX(curr)[,'FILTER'] == 'PASS'))/length(which(getFIX(curr)[,'FILTER'] == 'PASS'))
    },error=function(e){
	    print(e)
	    return(0)
    })
  }
  file.list <- clinical_file[,mutationFileField];
  get.file <- function(X) {
    tryCatch({
      na.exclude(X)[1]
    },error=function(e){NA})
  }
  file.list <- apply(file.list,1,get.file)
  mb <- unlist(lapply(file.list,getMutationCount,testdir=testdir))

  output <- data.frame(study=clinical_file$Study,
    patient=clinical_file$Patient,
    predictionscore=rescale(clinical_file$D_ISS,to=c(1,10))*rescale(clinical_file$D_Age,to=c(1,10))*rescale(mb,to=c(1,10)),
    highriskflag=clinical_file$D_ISS>2 & clinical_file$D_Age > 65 & mb > mutationburden)

  write.table(file=output.file,output,row.names=FALSE,sep="\t",quote=FALSE)
}

test.classifier()
