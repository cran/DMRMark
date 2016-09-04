dm_dmvnorm <- function(x,mean,sigma,log=FALSE,pd,logdet,invSigma) {
  if (missing(mean)) {
    mean <- matrix(0, ncol = ncol(x))
  }
  if(is.vector(mean)) {
    if(!is.null(dim(x))) {
      t_m = nrow(x)
    } else {
      t_m = 1
    }
  }
  if(missing(invSigma)) {
    if (missing(sigma)) {
      sigma <- diag(2)
    }
    invSigma <- tryCatch(chol2inv(chol(sigma)), 
                         error = function(e) {
                           print("sigma error") 
                           return(diag(2)/1000)
                         }
    )
  }
  if(missing(logdet)) {
    logdet <- log(det(sigma))
  }
  if(!missing(pd)) {
    y = reformData(x, pd)
    n = ncol(pd) - 1
  } else {
    y = x
    n = 1
  }
  
  dens <- rep(0, t_m)
  
  if (is.matrix(y)) {
    k <- ncol(y) - sum(is.na(colSums(y)))
  } else {
    k <- length(y) - sum(is.na(y))
    if (!is.na(y2[1])) {
      if (!is.na(y2[2])) {
        dens <- dmvnorm(y2, mean, sigma, TRUE)
      } else {
        dens <- dens + dnorm(y[1], mean[1], sqrt(sigma[1,1]), TRUE)
      }
    } else {
      dens <- dens + dnorm(y[2], mean[2], sqrt(sigma[2,2]), TRUE)
    }
    if (log) {
      return(dens)
    } else {
      return(exp(dens))
    }
  }

  
  if (!missing(pd)) {
    #get density mode
    for (nn in 1 : n) {
      y2 <- y[(nn*t_m-t_m+1):(nn*t_m),]
      #Only need to check the first row
      if (!is.na(y2[1,1])) {
        if (!is.na(y2[1,2])) {
          dens <- dens + dmvnorm(y2, mean, sigma, TRUE)
          
        } else {
          dens <- dens + dnorm(y2[,1], mean[1], sqrt(sigma[1,1]), TRUE)
        }
      } else {
        dens <- dens + dnorm(y2[,2], mean[2], sqrt(sigma[2,2]), TRUE)
      }
    }
    
  } else {
    # get likelihood mode
    i1 <- which(!is.na(y[,1]) & !is.na(y[,2]))
    i2 <- which(!is.na(y[,1]) & is.na(y[,2]))
    i3 <- which(is.na(y[,1]) & !is.na(y[,2]))
    
    dens[i1] <- dmvnorm(y[i1,], mean, sigma, TRUE)
    dens[i2] <- dnorm(y[i2,1], mean[1], sqrt(sigma[1,1]), TRUE)
    dens[i3] <- dnorm(y[i3,2], mean[2], sqrt(sigma[2,2]), TRUE)
  }
  
  if (log) {
    return(dens)
  } else {
    return(exp(dens))
  }
}

dnormP <- function(data, mu, Sigma, log = FALSE, pd) {
  n = ncol(pd) - 1
  y = reformData(data, pd)
  t_m = nrow(data)
  dens = rep(0, t_m)
  for (nn in 1 : n) {
    X <- y[(nn*t_m-t_m+1):(nn*t_m),]
    
    #Only need to check the first row
    if (!is.na(X[1,1])) {
      if (!is.na(X[1,2])) {
        dens <- dens + dnorm(X[,1], mu, sqrt(Sigma[1,1]), TRUE)
        dens <- dens + dnorm(X[,2]-X[,1], 0, sqrt(Sigma[2,2]), TRUE)
      } else {
        dens <- dens + dnorm(X[,1], mu, sqrt(Sigma[1,1]), TRUE)
      }
    } else {
      dens <- dens + dnorm(X[,2], mu, sqrt(Sigma[1,1]+Sigma[2,2]), TRUE)
    }
  }
  
  if (log) {
    return(dens)
  } else {
    return(exp(dens))
  }
}

getTrDens <- function(charL, L) {
  t_m <- length(L)
  k <- length(charL)
  trDens <- array(1/k,c(t_m,k,k))
  charL2 <- -1/charL
  rho <- 1/k*exp(matrix(L, ncol = 1) %*% charL2)
  for (i in 1 : k) {
    for (j in 1 : k) {
      if (j == 6) {
        trDens[,j,] <- rep(1/k,k)
        break
      }
      if (i == j) {
        trDens[,j,i] <- 1/k + (k-1)*rho[,j]
      } else {
        trDens[,j,i] <- 1/k - rho[,j]
      }
    }
  }
  trDens
}

getDens <- function(data, mu, sigma, pd) {
  k <- nrow(mu)
  t_m <- nrow(data)
  dens <- matrix(0,t_m,k)
  
  for (i in 1 : k) {
    if (i %in% c(1,2,5)) {      
      dens[,i] <- dnormP(data, mu[i,1], sigma[,,i], TRUE, pd)
    } else {
      dens[,i] <- dm_dmvnorm(data, mu[i,], sigma[,,i], TRUE, pd)
    }
  }
  dens
}

##########################
#Viterbi algorithm with log dens (private)

logViterbi <- function(init,A,B,ntimes) {
  # returns the most likely state sequence
  
  nt <- dim(B)[1]
  ns <- ncol(init)
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  const <- apply(B,1,max)  
  B <- B - const
  B <- exp(B)
  #B <- B / rowSums(B)
  ns <- dim(B)[2]  
  prior <- init
  
  delta <- psi <- matrix(nrow=nt,ncol=ns)
  state <- vector(length=nt)
  
  for(case in 1:lt) {
    # initialization
    delta[bt[case],] <- prior[case,]*B[bt[case],]
    delta[bt[case],] <- delta[bt[case],]/(sum(delta[bt[case],]))
    psi[bt[case],] <- 0
    # recursion
    if(ntimes[case]>1) {
      for(tt in ((bt[case]):(et[case]-1))) {
        for(j in 1:ns) {
          # A[tt,j,]: from tt to tt+1, from j to any
          # A[tt,,k]: from tt to tt+1, from any to k
          # delta[tt+1,j]: the prob of most likely path to j
          k <- which.max(delta[tt,]*A[tt,,j]) # if next is j, k is the most likely current state
          delta[tt+1,j] <- delta[tt,k]*A[tt,k,j]*B[tt+1,j]
          psi[tt+1,j] <- k
        }
        delta[tt+1,] <- delta[tt+1,]/(sum(delta[tt+1,]))
      }
    }
    
    # trace maximum likely state
    state[et[case]] <- which.max(delta[et[case],])
    if(ntimes[case]>1) {
      for(i in (et[case]-1):bt[case]) {
        state[i] <- psi[i+1,state[i+1]]
      }
    }      
  }
  
  return(list(states = state, delta = delta))
}
####################

###########################################################
# Private function for internal use
# Change to old version options for capatibility

oldOption <- function (opts) {
  k <- 4
  mu0_old <- matrix(0,k,2)
  mu0_old[1:2,1] <- opts$theta0
  mu0_old[3:4,] <- opts$mu0
  
  v0_old <- c(opts$alpha12N[1:2],opts$nu0)
  
  K0_old <- opts$kappa0
  
  A0_old <- array(0, dim = c(2,2,k))
  A0_old[1,1,1:2] <-opts$beta12N[1:2]
  A0_old[2,2,1] <- opts$beta12N[3]
  A0_old[2,2,2] <- opts$alpha12N[3]
  A0_old[,,3:4] <- opts$A0
  
  gap <- c(opts$D_mu[2], -opts$D_mu[1])
  
  old_opts <- list(alpha = opts$pi0,
                   cmu0 = opts$cmu0,
                   mu0 = mu0_old,
                   K0 = K0_old,
                   v0 = v0_old,
                   A0 = A0_old,
                   gap = gap,
                   lvl = (1-opts$chi_alpha),
                   burnin = opts$burnin,
                   sampling = (opts$nsamples*opts$sampleSep),
                   nsamples = opts$nsamples,
                   sampleSep = opts$sampleSep,
                   onHMM = opts$onHMM,
                   track = opts$track,
                   verbose = opts$verbose)
  
  return(old_opts)
}
###########################################################

###########################################################
# Private function for internal use
# Heuristics for faster computing

myHeuristic <- function (mv, pd, L, starting, opts) {
  
  init <- FullSample(mv, L, starting, opts$alpha, opts$mu0, opts$K0, opts$v0, opts$A0, pd, opts$cmu0, 
                     gap=opts$gap, lvl=opts$lvl, 
                     burnin = 10, sampling = 1, 
                     track = FALSE, withHMM = FALSE, 
                     verbose = FALSE)
  
  results <- FullSample(mv, L, starting, opts$alpha, opts$mu0, opts$K0, opts$v0, opts$A0, pd, opts$cmu0, 
                        gap=opts$gap, lvl=opts$lvl, 
                        burnin = 10, nsamples = 20, sampleSep = 3,
                        track = FALSE, withHMM = TRUE, 
                        verbose = FALSE, randomInit = FALSE)
  
  return(results)
}
############################################################

###########################################################
# Private function for internal use
# Change old-version parameters to new version

toNewPars <- function(pars) {
  
  if (length(dim(pars$mu)) == 2) {
    theta <- pars$mu[1:2,1]
    mu <- pars$mu[3:4,]
    
    sigma12 <- pars$sigma[1,1,1:2]
    sigmaN <- pars$sigma[2,2,1]
    Sigma34 <- pars$sigma[,,3:4]
    
    return(list(theta = theta,
                mu = mu,
                sigma12 = sigma12,
                sigmaN = sigmaN,
                Sigma34 = Sigma34,
                charL = pars$charL,
                init = pars$init))
  } else {
    theta <- pars$mu[1:2,1,]
    mu <- pars$mu[3:4,,]
    
    sigma12 <- pars$sigma[1,1,1:2,]
    sigmaN <- pars$sigma[2,2,1,]
    Sigma34 <- pars$sigma[,,3:4,]
    
    return(list(theta = theta,
                mu = mu,
                sigma12 = sigma12,
                sigmaN = sigmaN,
                Sigma34 = Sigma34,
                charL = pars$charL,
                init = pars$init))
  }
}

###########################################################