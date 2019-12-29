library(cvTools)
library(glmnet)
# Utilize UAMS dataset for training
training <- combine(eval.data$GSE24080UAMS[[1]]$training,eval.data$GSE24080UAMS[[1]]$test)
emc92.features <- rbind(result.data$emc92$GSE24080UAMS[[1]],test.result.data$emc92$GSE24080UAMS[[1]]);
uams70.features <- rbind(result.data$uams.70$GSE24080UAMS[[1]],test.result.data$uams.70$GSE24080UAMS[[1]]);
uams17.features <- rbind(result.data$uams.17$GSE24080UAMS[[1]],test.result.data$uams.17$GSE24080UAMS[[1]]);
uams5.features <- rbind(result.data$uams.5$GSE24080UAMS[[1]],test.result.data$uams.5$GSE24080UAMS[[1]]);
emc92.features <- emc92.features[match(sampleNames(training),emc92.features$ID),]
uams70.features <- uams70.features[match(sampleNames(training),uams70.features$ID),]
uams17.features <- uams17.features[match(sampleNames(training),uams17.features$ID),]
uams5.features <- uams5.features[match(sampleNames(training),uams5.features$ID),]
all.training.uams <- data.frame(PFStimeMonths=pData(training)[['PFStimeMonths']],
                                Progression=pData(training)[['Progression']],
                                emc92=emc92.features$raw.score,
                                uams70.features=uams70.features$raw.score,
                                uams17.features=uams17.features$raw.score,
                                uams5.features=uams5.features$raw.score)
all.training.uams <- all.training.uams[all.training.uams$PFStimeMonths > 0,]


folds <- cvFolds(n=nrow(all.training.uams),K=10)
# Function to minimize to identify the best-fit alpha value
train <- function(par,numeric.matrix,folds,classes,mc) {
  library(glmnet)
  library(doParallel)
  registerDoParallel(cl = mc)
  model <- cv.glmnet(y=classes,x=numeric.matrix,family = 'cox',foldid = folds$which,parallel = T)
  model$cvm[model$lambda == model$lambda.1se]
}
mc <- makeCluster(detectCores())
classes <- as(Surv(time=as.numeric(all.training.uams$PFStimeMonths),event=all.training.uams$Progression),'matrix')
numeric.matrix=as.matrix(all.training.uams[,-c(1,2)])
optimized <- optim(par=0,train,method='Brent',lower=0,upper=1,numeric.matrix=numeric.matrix,folds=folds,classes=classes,mc)
model <- cv.glmnet(x=numeric.matrix,y=classes,family='cox',alpha=optimized$par,foldid = folds$which,parallel = T)
stopCluster(mc)
coef <- predict(model,newx=numeric.matrix,s=model$lambda.1se,type='coefficients')