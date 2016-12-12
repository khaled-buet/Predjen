featurelist = colnames(features)
featurelist[length(featurelist)+1] = "Score"
formuli = featureformula(featurelist)
y = featurizeddata$Score
Fselection <- Boruta(x = features, y = y)