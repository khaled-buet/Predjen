findposition <- function(sequence, pattern, position) {
  unlist(lapply(sequence,function(s)
    if(position %in% gregexpr(pattern = pattern, toString(s))[[1]]){1}else{0}
  ))
}