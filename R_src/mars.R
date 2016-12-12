#' MARS Regression
#'
#' This function takes multivariate adaptive regression splines feature list, dataset, leave one out value and a value for n-fold cross-validation. Now, it outputs spearman correlation based on provided dataset.
#' @param featurelist a list of feature
#' @param featuredata provided dataset
#' @param leaveonegene 1 means this method will perform leave-one-out cross-validation; otherwise it will perform n-fold cross-validation
#' @param kfold a value for cross validation, by default it is set to 10
#' @return a vector of spearman correlations for all runs
#' @export
#' @examples
#' featurelist = c("X30mer", "Percent.Peptide", "Amino.Acid.Cut.position","predictions")
#' #suppose we have a file as '../crisprpred/data-raw/sample_data.csv' and current directory is set to '../crisprpred'
#' #setwd('..')
#' dir = getwd()
#' filepath = paste0(dir,'/data-raw/sample_data.csv')
#' data = read.csv(filepath)
#' mars(featurelist,data,0)
mars = function(featurelist, featuredata, kfold = 10) {
  fformula = featureformula(featurelist)
  modelE = earth(as.formula(fformula), featuredata, nfold = kfold)
  predictionsE = predict(modelE, featuredata)
  j = 1
  rmseE = c()
  spcor = c()
  step = length(featuredata[,1]) / kfold
  for (i in 1:(kfold - 1)) {
    errorE = featuredata$predictions[j:(j + step)] - predictionsE[j:(j + step)]
    rmseE = c(rmseE, rmse(errorE))
    spcor = c(spcor,cor(featuredata$predictions[j:(j + step)], predictionsE[j:(j +
                                                                                 step)], method = 'spearman'))
    j = j + step
  }
  cat("MARS Regression RMSE (n-fold):", rmseE, "\n")
  cat("Spearman Cor. for MARS (n-fold):", spcor, "\n")
  
  return(spcor)
}
