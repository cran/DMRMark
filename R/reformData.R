reformData <- function (mv, pd=NULL) {
  if (any(is.na(mv))) {
    #stop("Reform data error: NA not handled.")
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
  t_m = nrow(mv)
  mv2 <- matrix(nrow = 0, ncol = 2)
  for (i in 1 : n) {
    temp = which(pd[,i] == 1)
    if (length(temp) == 1) {
      if (pd[temp,n+1] == 0) {
        # Normal
        mv2 <- rbind(mv2, cbind(mv[,temp],rep(NA,t_m)))
        
      } else {
        # Tumor
        mv2 <- rbind(mv2, cbind(rep(NA,t_m),mv[,temp]))
        
      }
    } else {
      if (pd[temp[1],n+1] == 0) {
        # Normal first
        mv2 <- rbind(mv2, cbind(mv[,temp[1]],mv[,temp[2]]))
        
      } else {
        # Tumor first
        mv2 <- rbind(mv2, cbind(mv[,temp[2]],mv[,temp[1]]))
        
      }
    }
  }
  mv2
}