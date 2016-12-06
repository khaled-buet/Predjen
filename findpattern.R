countpattern <- function(sequence, pattern) {
  unlist(lapply(sequence,function(s)
    if(-1 %in% gregexpr(pattern = pattern, toString(s))[[1]]){0}else{length(
      gregexpr(pattern = pattern, toString(s))[[1]])}
    ))
}