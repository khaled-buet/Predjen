#' Illustration of Featurization
#'
#' This function takes dataset and a list of features as input and produce a features-wise dataset. The number of columns in returned dataset is equal to the number of features in featurelist.
#'
#' @param sequences provided as dataframe
#' @param string a list of aminoacids or nucleotides
#' @param seq sequence based features. by default it is true.
#' @param seqord highest number of sequence which will be considered together
#' @param pos position specific features. by default it is true.
#' @param posord highest number of sequence which will be considered together
#' @return a featurized dataframe
#' @export
#' @examples
#' input = list("ABCDEFGHABDAACBBDEBGGGHHH", "ABCBDBEBEBBBDBDBFDFDFGGHHEEFFEECCCD")
#' string = c("A", "BD")
#' featuredata = featurization(input, string, seq = TRUE, pos = FALSE)
#' featuredata
library(gtools)
featurization <-
  function(sequences, string, seq = TRUE, seqorder = 2, pos = TRUE, posorder = 7) {
    features = data.frame(1:length(sequences))
    colnames(features)[length(features)] = "Serial"
    if (seq == TRUE) {
      for (s in 1:seqorder) {
        permu = permutations(
          n = length(string), r = s, v = string, repeats.allowed = TRUE
        )
        for (i in 1:length(permu[,1])) {
          temp = countpattern(sequence = sequences, pattern = paste(permu[i,], collapse = ''))
          #cat(length(temp),permu[i,],"\n")
          features = data.frame(features, temp)
          colnames(features)[length(features)] = paste(permu[i,], collapse = '')
        }
      }
    }
    if (pos == TRUE)
    {
      minlength = min(unlist(lapply(sequences, function(s) {
        nchar(toString(s))
      })))
      for (p in 1:posorder) {
        permu = permutations(
          n = length(string), r = p, v = string, repeats.allowed = TRUE
        )
        for (i in 1:length(permu[,1])) {
          for (j in 1:minlength) {
            temp = findposition(sequence = sequences, pattern = paste(permu[i,], collapse = ''), j)
            cat(length(temp),permu[i,],"\n")
            features = data.frame(features, temp)
            colnames(features)[length(features)] = paste0(paste(permu[i,], collapse = ''), "_", j)
          }
        }
      }
    }
    features$Serial = NULL
    return(features)
  }