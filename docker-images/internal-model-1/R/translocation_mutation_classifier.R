library(vcfR)

args = commandArgs(trailingOnly=TRUE)

testdir = '/test-data'
outputdir = '/output'

test.classifier <- function(input.file=file.path(testdir,args[1]),
                            mutationburden=64.51,
                            mutationFileField=c('WES_mutationFileMutect','WES_mutationFileStrelkaSNV','RNASeq_mutationFileMutect'),
                            output.file=file.path(outputdir,'predictions.tsv')) {
  
  clinical_file <- read.csv(file=input.file,header=T,as.is=T,check.names=F)
  
  # Add a rudimentary imputation step
  clinical_file$D_ISS[is.na(clinical_file$D_ISS)] <- median(as.numeric(clinical_file$D_ISS),na.rm=T);
  clinical_file$CYTO_predicted_feature_08[is.na(clinical_file$CYTO_predicted_feature_08)] <- median(as.numeric(clinical_file$CYTO_predicted_feature_08),na.rm=T)
  clinical_file$CYTO_predicted_feature_05[is.na(clinical_file$CYTO_predicted_feature_05)] <- median(as.numeric(clinical_file$CYTO_predicted_feature_05),na.rm=T)
  if(all(is.na(clinical_file$D_ISS))) {
    clinical_file$D_ISS <- rep(3,length(clinical_file$D_ISS))
  }
  clinical_file$D_ISS <- as.numeric(clinical_file$D_ISS == 3)
  if(all(is.na(clinical_file$CYTO_predicted_feature_08))) {
    clinical_file$CYTO_predicted_feature_08 <- rep(0,length(clinical_file$D_ISS))
  }
  if(all(is.na(clinical_file$CYTO_predicted_feature_05))) {
    clinical_file$CYTO_predicted_feature_05 <- rep(0,length(clinical_file$D_ISS))
  }
  
  getMutationCount <- function(X,testdir) {
    tryCatch({
      curr <- read.vcfR(file.path(testdir,X),verbose=FALSE);
      length(which(grepl(pattern = '\\|TP53\\|',x = curr@fix[,8]) & grepl(pattern='missense',curr@fix[,8]) & getFIX(curr)[,'FILTER'] == 'PASS'))/length(which(grepl(pattern = '\\|TP53\\|',x = curr@fix[,8]) & getFIX(curr)[,'FILTER'] == 'PASS'))
    },error=function(e){0})
  }
  file.list <- clinical_file[,mutationFileField];
  get.file <- function(X) {
    tryCatch({
      na.exclude(X)[1]
    },error=function(e){NA})
  }
  file.list <- apply(file.list,1,get.file)
  mb <- unlist(lapply(file.list,getMutationCount,testdir=testdir))
  mb[is.na(mb)] <- median(na.exclude(mb))
  output <- data.frame(study=clinical_file$Study,
                       patient=clinical_file$Patient,
                       predictionscore=(mb+1)*(clinical_file$D_ISS+1)*(clinical_file$CYTO_predicted_feature_08+1)*(clinical_file$CYTO_predicted_feature_05+1),
                       highriskflag=((mb*100 > mutationburden & clinical_file$CYTO_predicted_feature_08 != 0) | (clinical_file$D_ISS != 0 & clinical_file$CYTO_predicted_feature_05 != 0))
  )
  
  write.table(file=output.file,output,row.names=FALSE,sep="\t",quote=FALSE)
}

test.classifier()
