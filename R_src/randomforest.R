#' Random Forest
#'
#' This function takes random forest feature list, dataset, leave one out value, a value for n-fold cross-validation, number of trees and sampline rate. Now, it outputs spearman correlation based on provided dataset.
#' @param featurelist a list of feature
#' @param featuredata provided dataset
#' @param leaveonegene 1 means this method will perform leave-one-out cross-validation; otherwise it will perform n-fold cross-validation
#' @param kfold a value for cross validation, by default it is set to 10
#' @param trees number of trees that are allowed to grow. default one is set to 500.
#' @param learningrate it is the sampling rate of data set. default is set to 0.6
#' @return a vector of spearman correlations for all runs
#' @export
#' @examples
#' featurelist = c("X30mer", "Percent.Peptide", "Amino.Acid.Cut.position","predictions")
#' #suppose we have a file as '../crisprpred/data-raw/sample_data.csv' and current directory is set to '../crisprpred'
#' #setwd('..')
#' dir = getwd()
#' filepath = paste0(dir,'/data-raw/sample_data.csv')
#' data = read.csv(filepath)
#' h2o.init()
#' randomforest(featurelist,data,leaveonegene=1)

randomforest = function(featurelist, featuredata, kfold = 10, trees = 50, learningrate = 0.6) {
  predict = featurelist[length(featurelist)]
  featurelist = featurelist[-length(featurelist)]
  preddata = as.vector(data[predict][,])
  h2o.init()
  inputdata.hex = as.h2o(featuredata, destination_frame = "inputdata.hex")
  randmodel = h2o.randomForest(
    x = featurelist, y = predict, training_frame = inputdata.hex, nfolds = kfold, ntrees = trees, sample_rate = learningrate
  )
  rfpred = h2o.predict(randmodel, inputdata.hex)
  #cat('DP:',as.vector(dppred),'\n')
  j = 1
  rmseE = c()
  spcor = c()
  step = length(featuredata[,1]) / kfold
  for (i in 1:(kfold - 1)) {
    rmseE = c(rmseE,rmse(preddata[j:(j + step)] - as.vector(rfpred)[j:(j + step)]))
    spcor = c(spcor,cor(preddata[j:(j + step)], as.vector(rfpred)[j:(j + step)], method = 'spearman'))
    j = j + step;
  }
  cat("Random Forest RMSE (n-fold):", rmseE,"\n")
  cat("Spearman Cor. for RF (n-fold):", spcor, "\n")
  return(spcor)
}
