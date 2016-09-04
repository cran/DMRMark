mvScatter <- function(mv, isDML, pd=NULL, nPlot=5000) {
  plot(-6.5:6.5,-6.5:6.5, type = "n",
       xlim = c(-6.5, 6.5), xlab = colnames(mv)[1],
       ylim = c(-6.5, 6.5), ylab = colnames(mv)[2])
  pch <- c(1,3)
  clr <- c("green","red")
  lwd <- c(1,1)
  
  t_m2 <- nrow(mv)
  if (t_m2 > nPlot) {
    temp <- sample(1:t_m2, nPlot)
    trueSig <- isDML[temp]
    mv <- mv[temp,]
    t_m2 <- nPlot
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
  
  prb <- list()
  prb[[1]] <- which(isDML == 0)
  prb[[2]] <- which(isDML == 1)
  
  mv2 <- reformData(mv, pd)
  
  for (tp in 1 : length(pch)) {
    points(mv2[prb[[tp]],1],mv2[prb[[tp]],2], pch = pch[tp], col = clr[tp], lwd = lwd[tp])
  }
  lines(x = -7:7, y = -7:7, col = "blue", lwd = 2)
}