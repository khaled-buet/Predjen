#amins = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
library(e1071)
library(ROCR)
setwd('..')
dir = getwd()
print(dir)
source(paste0(dir,"/R_src/datapartition/datapartition.R"))
source(paste0(dir,"/R_src/featureformula/featureformula.R"))
setwd(paste0(dir,"/data-raw/"))
randomforestfeatures = read.csv('randomForestFeatures.csv')
fcres = read.csv('FC_plus_RES.csv')
features = c()
for(i in 1:length(randomforestfeatures[,2])){
  if(randomforestfeatures[i,2] > 0.1758){
    features = c(features,toString(randomforestfeatures[i,1]))
  }
}
features[length(features) + 1] = "geneThreshold"
cat('Reading features complete: total',  length(features), ' features ...\n')
ontargetdata = read.csv('onTargetFinalData.csv')
columns = names(ontargetdata)
for(i in 1:length(columns)){
  if(!columns[i] %in% features){
    ontargetdata[columns[i]] = NULL
  }
}
gene = fcres$Target.gene
ontargetdata = data.frame(ontargetdata, gene)
cat('Reading featured data complete ...\n')

genes = c(
  'CCDC101', 'CD13','CD15', 'CD28' ,'CD33' ,'CD43' ,'CD45' ,'CD5' ,'CUL3', 'H2-K', 'HPRT1' ,'MED12', 'NF1' ,'NF2', 'TADA1', 'TADA2B', 'THY1'
)
data = ontargetdata
predValue = c()
for (i in 1:length(genes)) {
  dataset = datapartition(
    data, targetcolumn = "gene", leaveonegene = TRUE, genename = genes[i]
  )
  training = dataset[[1]]
  testing = dataset[[2]]
  bformula = featureformula(features)
  svmmodel = svm(
    as.formula(bformula), training
  )
  svmpred = predict(svmmodel, testing)
  svmpred = as.vector(svmpred)
  predValue = c(predValue, svmpred)
  cat("Running SVM iteration no. = ",i,"  ...\n")
}
