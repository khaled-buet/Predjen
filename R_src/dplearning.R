#' Deep Learning
#'
#' This function takes deep learning feature list, dataset, leave one out value, a value for n-fold cross-validation and learning rate. Now, it outputs spearman correlation based on provided dataset.
#' @param featurelist a list of feature
#' @param featuredata provided dataset
#' @param leaveonegene 1 means this method will perform leave-one-out cross-validation; otherwise it will perform n-fold cross-validation
#' @param kfold a value for cross validation, by default it is set to 10
#' @param learningrate learning rate value. Default one is set to 0.6
#' @return a vector of spearman correlations for all runs
#' @export
#' @examples
#' featurelist = c("X30mer", "Percent.Peptide", "Amino.Acid.Cut.position","predictions")
#' #suppose we have a file as '../crisprpred/data-raw/sample_data.csv' and current directory is set to '../crisprpred'()
#' setwd('..')
#' dir = getwd()
#' filepath = paste0(dir,'/data-raw/sample_data.csv')
#' data = read.csv(filepath)
#' h2o.init()
#' dplearning(featurelist, data,leaveonegene=1)
dplearning = function(featurelist,featuredata, kfold = 10, learningrate = 0.6) {
  predict = featurelist[length(featurelist)]
  featurelist = featurelist[-length(featurelist)]
  preddata = as.vector(data[predict][,])
  h2o.init()
  inputdata.hex = as.h2o(featuredata, destination_frame = "inputdata.hex")
  dpmodel = h2o.deeplearning(
    x = featurelist, y = predict, training_frame = inputdata.hex, nfolds = kfold, rate = learningrate
  )
  dppred = h2o.predict(dpmodel, inputdata.hex)
  #cat('DP:',as.vector(dppred),'\n')
  j = 1
  rmseE = c()
  spcor = c()
  step = length(featuredata[,1]) / kfold
  for (i in 1:(kfold - 1)) {
    rmseE = c(rmseE,rmse(preddata[j:(j + step)] - as.vector(dppred)[j:(j + step)]))
    spcor = c(spcor,cor(preddata[j:(j + step)], as.vector(dppred)[j:(j + step)], method = 'spearman'))
    j = j + step;
  }
  cat("Deep Learning RMSE (n-fold):", rmseE,"\n")
  cat("Spearman Cor. for DL (n-fold):", spcor, "\n")
  return(spcor)
}
