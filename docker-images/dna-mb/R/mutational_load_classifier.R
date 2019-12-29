library(vcfR)

args = commandArgs(trailingOnly=TRUE)

testdir = '/test-data'
outputdir = '/output'

test.classifier <- function(input.file=file.path(testdir,args[1]),
mutationburden=64.51,
mutationFileField=c('WES_mutationFileMutect','WES_mutationFileStrelkaSNV','RNASeq_mutationFileMutect'),output.file=file.path(outputdir,'predictions.tsv')) {
  clinical_file <- read.csv(file=input.file,header=T,as.is=T,check.names=F)

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
    predictionscore=mb,
    highriskflag=mb > mutationburden)

  write.table(file=output.file,output,row.names=FALSE,sep="\t",quote=FALSE)
}

test.classifier()
