read.data = function(){
  d = read.table("exogenic.txt", header = TRUE, fill = TRUE)
  d = head(d, nrow(d) - 3)
  for(i in 1:3){
    d[, i] = as.numeric(d[, i])
  }
  return(d)
}

