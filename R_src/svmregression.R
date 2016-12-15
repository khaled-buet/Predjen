#' SMV Regression
#'
#' This function takes svm regression formula, dataset, leave one out value  and a value for n-fold cross-validation. Now, it outputs spearman correlation based on provided dataset.
#' @param featurelist a list of feature
#' @param featuredata provided dataset
#' @param leaveonegene 1 means this method will perform leave-one-out cross-validation; otherwise it will perform n-fold cross-validation
#' @param kfold a value for cross validation, by default it is set to 10
#' @return a vector of spearman correlations for all runs
#' @export
#' @examples
#' featurelist = c("X30mer", "Percent.Peptide", "Amino.Acid.Cut.position","predictions")
#' #suppose we have a file as '../crisprpred/data-raw/sample_data.csv' and current directory is set to '../crisprpred'
#' #set('..')
#' dir = getwd()
#' filepath = paste0(dir,'/data-raw/sample_data.csv')
#' data = read.csv(filepath)
#' svmregression(featurelist,data,1)
svmregression = function(featurelist,featuredata, kfold = 10) {
  fformula = featureformula(featurelist)
  modelS = svm(as.formula(fformula), featuredata, cross = kfold)
  predictionsS = predict(modelS, featuredata)
  j = 1
  rmseS = c()
  spcor = c()
  step = length(featuredata[,1]) / kfold
  for (i in 1:(kfold - 1)) {
    errorS = featuredata$Score[j:(j + step)] - predictionsS[j:(j + step)]
    rmseS = c(rmseS, rmse(errorS))
    spcor = c(spcor, cor(featuredata$Score[j:(j + step)], predictionsS[j:(j +
                                                                                  step)], method = 'spearman'))
    cat(i,": SVM SPCOR = ", spcor, "\n")
    j = j + step
  }
  cat("SVM Regression RMSE (n-fold):", rmseS, "\n")
  cat("Spearman Cor. for SVM (n-fold):", spcor, "\n")
  
  return(spcor)
}
