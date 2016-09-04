reformData <- function (mv, pd=NULL) {
  if (any(is.na(mv))) {
    stop("Reform data error: NA not handled.")
  }
  
  if (is.null(pd)) {
    n <- ncol(mv)/2
    if (n > 1) {
      patient <- factor(c(1:n,1:n))
      type = c(rep("Normal",n),rep("Tumour",n))
      pd <- model.matrix(~patient + type + 0)
    } else {
      pd <- matrix(0,2,2)
      pd[1:2,1] <- 1
      pd[2,2] <- 1
    }
  }
  
  n <- ncol(pd) - 1
  t <- ncol(mv)
  if (is.null(t)) {
    t <- length(mv)
  }
  
  mv2 <- NULL
  for (i in 1 : n) {
    b <- matrix(0, t, 2)
    temp = which(pd[,i] == 1)
    if (length(temp) == 1) {
      if (pd[temp,n+1] == 0) {
        b[temp, 1] <- 1
        b[temp, 2] <- NA
      } else {
        b[temp, 2] <- 1
        b[temp, 1] <- NA
      }
    } else {
      b[temp[1], pd[temp[1],n+1]+1] <- 1
      b[temp[2], pd[temp[2],n+1]+1] <- 1
    }
    mv2 <- rbind(mv2, mv %*% b)
  }
  mv2
}