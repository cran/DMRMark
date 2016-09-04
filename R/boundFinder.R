boundFinder <- function(mv, prop=0.1) {
  n <- ncol(mv)/2
  if (n > 1) {
    temp <- cbind(rowMeans(mv[,1:n]), rowMeans(mv[,(1:n)+n]))
  } else {
    temp <- mv
  }
  temp <- temp[,1] - temp[,2]
  temp <- quantile(temp,(1:1000)*0.001, na.rm = TRUE)
  i <- 1
  j <- which(temp > -temp[i])
  if (length(j) == 0) {
    j <- 1000
  } else {
    j <- min(which(temp > -temp[i]))
  }
  
  while ((i+1000-j) < (prop*1000)) {
    i <- i + 1
    j <- which(temp > -temp[i])
    if (length(j) == 0) {
      j <- 1000
    } else {
      j <- min(which(temp > -temp[i]))
    }
  }
  gap <- c(temp[i], temp[j])
  gap
}