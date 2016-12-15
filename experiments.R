featurelist = colnames(features)
featurelist[length(featurelist) + 1] = "Score"
formuli = featureformula(featurelist)
y = featurizeddata$Score
Fselection <- Boruta(x = features, y = y)
importantFeature = unlist(lapply(Fselection$finalDecision, function(a)
  if (a == 'Confirmed' || a == 'Tentative') {
    a
  }))
tentativeFeature = unlist(lapply(Fselection$finalDecision, function(a)
  if (a == 'Tentative') {
    a
  }))

p = ap$`Pr(>F)`
f = row.names(ap)
p = p[-length(p)]
f = f[-length(f)]
relevantfeatures = c()
for (i in 1:length(p)) {
    if (as.numeric(p[i]) < 0.0001) {
      cat(i,":",p[i],"\n")
      relevantfeatures[length(relevantfeatures) + 1] = f[i]
    }
}
f = relevantfeatures
f[length(f) + 1] = "Score"
l = lmregression(f, featurizeddata)
s = svmregression(f, featurizeddata)
r = randomforest0(f, featurizeddata)
dd = data.frame(l,s,r)
colnames(dd) = c('LR','SVM','RF')
matplot(dd,type = c("l","b","o"))