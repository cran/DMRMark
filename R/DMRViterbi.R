DMRViterbi <- function(mv, pars, L=rep(1,nrow(mv)), starting=NULL, pd=NULL, 
                       region=TRUE, orderBy=c("max", "mean", "min"), VitP=NULL) {
  
  if (is.null(starting)) {
    temp <- which((L > 100000) | L < 0)
    # The starting positions of new chains
    starting <- c(1, temp[-length(temp)]+1)
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
  
  t_m <- length(L)
  ntimes <- diff(c(starting, t_m+1))
  
  ###############################
  # Change to old-version parameters
  mu <- matrix(0,4,2)
  mu[1:2,1] <- pars$theta
  mu[3:4,] <- pars$mu
  
  Sigma <- array(0, dim = c(2,2,4))
  Sigma[1,1,1:2] <- pars$sigma12
  Sigma[2,2,1:2] <- pars$sigmaN
  Sigma[,,3:4] <- pars$Sigma34
  ###############################
  
  dens <- getDens(mv, mu, Sigma, pd)
  trDens <- getTrDens(pars$charL, L)
  init <- matrix(rep(pars$init, t_m), t_m, byrow = TRUE)
  
  fm2 <- logViterbi(init, trDens, dens, ntimes)
  if (!region) {
    return(fm2)
  }
  
  states <- apply(fm2$delta, 1, which.max)
  regBegin <- rep(FALSE,t_m)
  regBegin[c(1,t_m)] <- TRUE
  regBegin[starting] <- TRUE
  regBegin[which(diff(states) != 0)+1] <- TRUE
  regBound = which(regBegin)
  regNum <- length(regBound) - 1
  
  results <- data.frame(begin=regBound[-regNum],
                        ending=regBound[-1]-1,
                        MAP_state=rep(0,regNum),
                        minVP=rep(-1,regNum),
                        meanVP=rep(-1,regNum),
                        maxVP=rep(-1,regNum))

  for (i in 1 : regNum) {
    b <- results$begin[i]
    e <- results$ending[i]
    st <- states[regBound[i]]
    
    results$MAP_state[i] <- st
    results$minVP[i] <- min(fm2$delta[b:e,st])
    results$meanVP[i] <- exp(mean(log(fm2$delta[b:e,st]), na.rm = TRUE))
    results$maxVP[i] <- max(fm2$delta[b:e,st])
  }

  #order by VP decreasingly, from DM to non-DM
  VP <- results$maxVP
  temp <- match.arg(orderBy)
  if (temp == "mean") {
  VP <- results$meanVP
  } else if (temp == "min") {
    VP <- results$minVP
  }

  results <- results[order(VP,  decreasing = TRUE),]
  isNDM <- (results$MAP_state <= 2) + 0
  results <- results[order(isNDM),]

  return(results)
}
