library(ROCR)
drawROC = function(TrueObservation, PredictedValue)
{
  pred <- prediction(PredictedValue, TrueObservation);
  ROC.perf <- performance(pred, "tpr", "fpr");
  plot (ROC.perf);
  
}