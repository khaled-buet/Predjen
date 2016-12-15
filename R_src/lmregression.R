#' Linear Regression
#'
#' This function takes featurelist, dataset, leave-one-out value  and a value for n-fold cross-validation. It creates a formula and outputs spearman correlation based on provided dataset.
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
#' lmregression(featurelist,data,1)
lmregression = function(featurelist, featuredata, kfold = 10) {
  fformula = featureformula(featurelist)
  model1 = lm(as.formula(fformula), featuredata)
  predictionsS = stats::predict(model1, featuredata)
  j = 1
  rmseS = c()
  spcor = c()
  step = length(featuredata[,1]) / kfold
  for (i in 1:(kfold - 1)) {
    errorS = featuredata$Score[j:(j + step)] - predictionsS[j:(j + step)]
    rmseS = c(rmseS, rmse(errorS))
    spcor = c(spcor, stats::cor(featuredata$Score[j:(j + step)], predictionsS[j:(j +
                                                                                         step)], method = 'spearman'))
    cat(i,": LM SPCOR = ", spcor, "\n")
    j = j + step
  }
  cat("LR Regression RMSE (n-fold):", rmseS, "\n")
  cat("Spearman Cor. for LR (n-fold):", spcor, "\n")
  return(spcor)
}
