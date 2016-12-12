data = read.csv('sequencefeaturized.csv')
aminoacid = c(
  "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"
)
data = data.frame(data)
features = featurization(sequences = data$Sequence, string = aminoacid)
featurizeddata = data.frame(data, features)
write.csv(featurizeddata, 'finaldata.csv')

