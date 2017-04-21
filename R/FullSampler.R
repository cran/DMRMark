FullSample <- function(data, L, starting, 
                       alpha, mu0, K0, v0, A0, pd, cmu0, 
                       gap = c(0,0), lvl = 0.8,
                       burnin = 500, sampling = NULL, 
                       track = FALSE, withHMM = TRUE, 
                       verbose = FALSE, nsamples = 100, sampleSep = 10,
                       randomInit = TRUE) {
  
  k <- length(alpha)
  t_m <- nrow(data)
  states <- sample(1:k, t_m, replace = TRUE)
  charL0 <- rep(100000,k) # 
  for (kk in 1 : k) {
    charL0[kk] <- exp(rnorm(1, cmu0[kk], 0.1))
  }
  if (!withHMM) {
    charL0 <- rep(0.1,k)
  }
  
  if(any(is.na(data))) {
    stop("NA not handled in input data.")
  }
  
  temp <- reformData(data, pd)
  a1 <- any(is.na(temp[,1]))
  a2 <- any(is.na(temp[,2]))
  if (a1 & a2) {
    cat("Unpaired samples detected in both groups, estimation can be slow.\n")
    noPair <- TRUE
  } else {
    noPair <- FALSE
  }
  rm(temp)
  
  if(!is.null(sampling)) {
    # mask following two parameters
    sampleSep <- 10
    nsamples <- ceiling(sampling/sampleSep)
  } else {
    sampling <- sampleSep * nsamples
  }
  
  # For storing results
  states_map <- matrix(0, t_m, k)
  if (track) {
    mu_map <- array(0,c(k,2,burnin+sampling+1))
    sigma_map <- array(0,c(2,2,k,burnin+sampling+1))
    init_map <- matrix(0, burnin+sampling+1,k)
    charL_map <- matrix(0,burnin+sampling+1,k)
  } else {
    mu_map <- matrix(0,k,2)
    sigma_map <- array(0,c(2,2,k))
    init_map <- rep(0, k)
    charL_map <- rep(0,k)
    coef <- 0
  }
  
  mu <- matrix(0,k,2)
  sigma <- array(0, c(2,2,k))
  if (randomInit) {
    for (kk in 1 : k) {
      mu[kk,] <- rnorm(2,mu0[kk,],1/K0[kk])
    }
    
    for (kk in 1:k) {
      if (kk == 3 | kk == 4) {
        sigma[,,kk] <- riwish(v0[kk], A0[,,kk])
      } else {
        sigma[1,1,kk] <- rinvgamma(1, v0[kk], A0[1,1,kk]) 
        if (kk == 1) {
          sigma[2,2,kk] <- rinvgamma(1, A0[2,2,2], A0[2,2,1]) 
        } else {
          sigma[2,2,kk] <- sigma[2,2,1]
        }
      }
    }
  } else {
    mu <- mu0
    for (kk in c(1,2)) {
      sigma[1,1,kk] <- A0[1,1,kk]
    }
    sigma[2,2,1:2] <- A0[2,2,1]
    
    for (kk in c(3,4)) {
      sigma[,,kk] <- A0[,,kk]
    }
  }

  initDens <- initSampler(states, starting, alpha)
  if (withHMM) {
    trDens <- getTrDens(charL0, L)
  }
  logDens <- getDens(data, mu0, sigma, pd)
  if (withHMM) {
    states <- HMMSampler(logDens, trDens, initDens, starting) #
  } else {
    states <- Ssampler(logDens, initDens)# 
  }

  for (b in 1 : burnin) {

    if (b %% 10 == 1) {
      #cat("burn-in at itr. ", b, "\n")
      if (verbose) {
        cat("burn-in at itr. ", b, "\n")
        print(mu)
      }
      #print(sigma[,,3])
    }
    
    initDens <- initSampler(states, starting, alpha)
    if (noPair) {
      parameter <- NIWSamplerP(data, states, mu0, K0, v0, A0, pd, gap, mu, sigma, lvl)
    } else {
      parameter <- NIWSampler(data, states, mu0, K0, v0, A0, pd, gap, mu, sigma, lvl)
    }

    mu <- parameter$mu
    sigma <- parameter$sigma
    if (withHMM) {
      trResults <- charLSampler(states, L, starting, cmu0, charL0)
      charL0 <- trResults$charL
      trDens <- getTrDens(charL0, L)
    }
    logDens <- getDens(data, mu, sigma, pd)
    if (withHMM) {
      states <- HMMSampler(logDens, trDens, initDens, starting) # 
    } else {
      states <- Ssampler(logDens, initDens)# 
    }

    if (track) {
      mu_map[,,b] <- mu
      sigma_map[,,,b] <- sigma
      init_map[b,] <- initDens
      charL_map[b,] <- charL0
    }
  }
  
  for (b in 1 : (sampling+1)) {
    if (b %% sampleSep == 1) {
      #cat("sampling at itr. ", b, "\n")
      if (verbose) {
        cat("sampling at itr. ", b, "\n")
        print(mu)
      }
      #print(sigma[,,3])
      for (tt in 1 : t_m) {
        states_map[tt,states[tt]] <- states_map[tt,states[tt]] + 1
      }
      if (!track) {
        mu_map <- mu_map + mu
        sigma_map <- sigma_map + sigma
        init_map <- init_map + initDens
        charL_map <- charL_map + log(charL0)
        coef <- coef + 1
      }
    }
    if (withHMM) {
      initDens <- initSampler(states, starting, alpha)
    } else {
      initDens <- initSampler(states, 1:t_m, alpha)
    }
    
    if (noPair) {
      parameter <- NIWSamplerP(data, states, mu0, K0, v0, A0, pd, gap, mu, sigma, lvl)
    } else {
      parameter <- NIWSampler(data, states, mu0, K0, v0, A0, pd, gap, mu, sigma, lvl)
    }
    
    mu <- parameter$mu
    sigma <- parameter$sigma
    if (withHMM) {
      trResults <- charLSampler(states, L, starting, cmu0, charL0)
      charL0 <- trResults$charL
      trDens <- getTrDens(charL0, L)
    }
    logDens <- getDens(data, mu, sigma, pd)
    if (withHMM) {
      states <- HMMSampler(logDens, trDens, initDens, starting) # 
    } else {
      states <- Ssampler(logDens, initDens)# 
    }
    
    if (track) {
      mu_map[,,burnin+b] <- mu
      sigma_map[,,,burnin+b] <- sigma
      init_map[burnin+b,] <- initDens
      charL_map[burnin+b,] <- charL0
    }
  }
  if (!track) {
    mu_map <- mu_map/coef
    sigma_map <- sigma_map/coef
    init_map <- init_map/coef
    charL_map <- exp(charL_map/coef)
  }
  
  return(list(states=states_map, mu=mu_map, sigma=sigma_map, init=init_map, charL = charL_map))
}

charLSampler <- function(states, L, starting, cmu0, charL0) {
  t_m <- length(states)
  skip <- c(starting[-1]-1,t_m)
  
  #######
  k <- length(charL0)
  Lt <- L[-skip]
  
  transLL <- function (At, Lt, cL) {
    rho <- 0.25*exp(-Lt/cL)
    tr <- cbind(0.25 + rho*3, 0.75 - rho*3)
    LogLik <- sum(log(tr) * At)
    
    LogLik
  }
  
  charL_all <- rep(0,k)
  LogLik_all <- rep(0,k)
  accept <- rep(FALSE, k)
  for (i in 1 : k) {
    if (i == 6) {
      charL_all[i] <- 0.10101
      next
    }
    
    At <- (states[-t_m] == i & states[-1] == i) + 0
    At0 <- (states[-t_m] == i & states[-1] != i) + 0
    At <- At[-skip]
    At0 <- At0[-skip]
    At <- cbind(At, At0)
    if (sum(At + At0) > 0) {
      logCL <- rnorm(1, mean = log(charL0[i]), sd = 0.15)
      charL <- exp(logCL)
#       print(transLL(At, Lt, charL))
#       print(dnorm(logCL, cmu0[i], 10/length(Lt), TRUE))
#       print(transLL(At, Lt, charL0[i]))
#       print(dnorm(log(charL0[i]), cmu0[i], 100/length(Lt), FALSE))
      LogLik <- transLL(At, Lt, charL) + dnorm(logCL, cmu0[i], 10, TRUE)
      LLold <- transLL(At, Lt, charL0[i]) + dnorm(log(charL0[i]), cmu0[i], 10, TRUE)
      acc <- exp(LogLik - LLold)
    } else {
      LLold <- -987654321
      acc <- -1
    }
    if (runif(1) < acc) {
      # Accept
      charL_all[i] <- charL
      LogLik_all[i] <- LogLik
      accept[i] <- TRUE
      
    } else {
      # rejest
      charL_all[i] <- charL0[i]
      LogLik_all[i] <- LLold
      accept[i] <- FALSE
    }
  }
  return(list(charL=charL_all, LogLik=LogLik_all, accept=accept))
}

initSampler <- function(states, starting, alpha) {
  k <- length(alpha)
  temp <- states[starting]
  n_k <- rep(0, k)
  for (kk in 1 : k) {
    n_k[kk] <- sum(temp == kk)
  }
  initDens <- rdirichlet(1, n_k+alpha)
  
  initDens
}

NIWSampler <- function(data, states, mu0, K0, v0, A0, pd, gap, mu_old = NULL, sigma_old = NULL, lvl = 0.8) {
  k <- nrow(mu0)
  mu_all <- matrix(0, k, 2)
  sigma_all <- array(rep(c(1,0,0,1),k), dim = c(2,2,k))
  for (kk in c(1,2)) {
    # For empty cluster
    if(sum(states == kk) <= 1) {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
      mu_all[kk,] <- mu_old[kk,]# rnorm(1, mu0[kk,], sigma_all[,,kk]/K0[kk])
      # cat("strange: ", kk, "\n")
      next
    }  
    
    # First compute the sigma of (C-N)
    # only for 1,2,5
    if (kk == 1) {
      X <- data[states %in% c(1,2,5),]
      X <- reformData(X, pd)
      X <- X[!is.na(rowSums(X)),]
      n <- nrow(X)
      X <- X[,2] - X[,1]
      S <- sum(X^2)
      sigma_all[2,2,kk] <- rinvgamma(1, A0[2,2,2]+n/2, A0[2,2,1]+S/2)
      
    } else if (kk == 2 | kk == 5) {
      sigma_all[2,2,kk] <- sigma_all[2,2,1]
    }
    
    # Compute the sigma & mu of N or C
    # only for 1,2,5
    X <- data[states == kk,]
    X <- reformData(X, pd)
    
    a1 <- any(is.na(X[,1]))
    a2 <- any(is.na(X[,2]))
    if (!a1) {
      X <- X[,1]
    } else if (!a2){
      X <- X[,2]
    } else {
      stop("Unexpected NA produced, please check pd.")
    }
    
    mu <- mean(X)
    if (kk == 5) {
      mu <- 0
      mu_all[kk,1] <- 0
      S <- sum(X^2)/2
      sigma_all[1,1,kk] <- rinvgamma(1, v0[kk]+n/2, A0[1,1,kk]+S)
    } else {
      n <- length(X)
      Kn <- K0[kk] + n
      S <- (n*K0[kk])/Kn*sum((mu - mu0[kk,1])^2)/2 + sum((X - mu)^2)/2
      sigma_all[1,1,kk] <- rinvgamma(1, v0[kk]+n/2, A0[1,1,kk]+S)
      mu_n <- (K0[kk]*mu0[kk,1] + n*mu) / Kn
      mu_all[kk,1] <- rnorm(1, mu_n, sqrt(sigma_all[1,1,kk]/Kn))
    }
  } 
    
  for (kk in c(3,4)) {
    # For empty cluster
    if(sum(states == kk) <= 1) {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
      mu_all[kk,] <- mu_old[kk,]# rnorm(1, mu0[kk,], sigma_all[,,kk]/K0[kk])
      # cat("strange: ", kk, "\n")
      next
    } 
    
    # compute the SIGMA of DML
    # only for 3,4
    X <- data[states == kk,]
    X <- reformData(X, pd)
    
    if (is.null(mu_old)) {
      i1 <- which(!is.na(X[,1]) & !is.na(X[,2]))
      n <- length(i1)
      mu <- colMeans(X, na.rm = TRUE)
      S <- var(X[i1,]) * (n-1)
      Kn <- K0[kk] + n
      vn <- v0[kk] + n
      mu_n <- (K0[kk] * mu0[kk,] + n * mu) / Kn
      mu_all[kk,] <- rmvnorm(1, mu_n, sigma_old[,,kk]/Kn)
      
    } else {
      # Truncated Normal Prior
      # Use MH sampling      
      
      acc <- -1
      mu_all[kk,] <- rmvnorm(1, mu_old[kk,], diag(2)*0.1)

      oldLL <- sum(dm_dmvnorm(X, mu_old[kk,], sigma_old[,,kk], TRUE)) + 
        dmvnorm(mu_old[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
      
      if (kk == 4) {
        #Hypo DMR
        if ((mu_all[kk,1] - mu_all[kk,2]) > gap[1]) {
          newLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) +
            dmvnorm(mu_all[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
          acc <- exp(newLL - oldLL)
        } else {
          if ((mu_old[kk,1] - mu_old[kk,2]) < gap[1]) {
            acc <- 0.5
          } else {
            acc <- -1
          }
        }
      } else if (kk == 3){
        #Hyper DMR
        if ((mu_all[kk,1] - mu_all[kk,2]) < -gap[2]) {
          newLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) +
            dmvnorm(mu_all[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
          acc <- exp(newLL - oldLL)
        } else {
          if ((mu_old[kk,1] - mu_old[kk,2]) > -gap[2]) {
            acc <- 0.5
          } else {
            acc <- -1
          }
        }
      } 
      if (acc < runif(1)) {
        mu_all[kk,] <- mu_old[kk,]
      }
    }
    
    oldLL <- sum(dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) + 
      log(diwish(sigma_old[,,kk], v0[kk], A0[,,kk]))
    
    s1_old <- sqrt(sigma_old[1,1,kk])
    s2_old <- sqrt(sigma_old[2,2,kk])
    #rho_old <- tan(sigma_old[1,2,kk]/abs(s1_old*s2_old)*pi/2)
    rho_old <- sigma_old[1,2,kk]/abs(s1_old*s2_old)
    
    s1 <- rnorm(1, s1_old, 0.025)
    s2 <- rnorm(1, s2_old, 0.025)
    rho <- rnorm(1, rho_old, 0.01)
    while (rho <= -1 | rho >= 1) {
      rho <- rnorm(1, rho_old, 0.025)
    }
    sigma_all[1,1,kk] <- s1^2
    sigma_all[2,2,kk] <- s2^2
    #sigma_all[1,2,kk] <- abs(s1*s2)*tanh(rho)*2/pi
    sigma_all[1,2,kk] <- abs(s1*s2)*rho
    sigma_all[2,1,kk] <- sigma_all[1,2,kk]

    pts <- ellipse(sigma_all[,,kk], centre = mu_all[kk,], level = lvl)
    dd <- pts[,2] - pts[,1]
    # New cov is ok
    if (all(dd > 0) | all(dd < 0)) {
      newLL <- sum(dmvnorm(X, mu_all[kk,], sigma_all[,,kk], TRUE)) + 
        log(diwish(sigma_all[,,kk], v0[kk], A0[,,kk]))
      acc <- exp(newLL - oldLL)
      if (acc >= runif(1)) {
        next
      }
    }

    # New cov is not ok, check if old cov is ok
    pts <- ellipse(sigma_old[,,kk], centre = mu_all[kk,], level = lvl)
    dd <- pts[,2] - pts[,1]
    if (all(dd > 0) | all(dd < 0)) {
      sigma_all[,,kk] <- sigma_old[,,kk]
    } else {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
    }
  }

  return(list(mu = mu_all, sigma = sigma_all))
}

HMMSampler <- function(logDens, trDens, initDens, starting = 1) {
  t_m <- nrow(logDens)
  lt <- length(starting)
  k = ncol(logDens)
  ending <- c(starting[-1]-1,t_m)
  
  const <- apply(logDens,1,max)
  logDens <- logDens - const
  dens <- exp(logDens)
  
  states <- rep(0, t_m)
  beta <- matrix(nrow=t_m,ncol=k)
  beta[ending,] <- dens[ending,]
  beta[ending,] <- beta[ending,] / rep(rowSums(beta[ending,]),k)
  
  for (case in 1 : lt) {
    if (ending[case] == starting[case]) {
      # chain with only 1 states
      next
    }
    for (i in (ending[case]-1) : starting[case]) {
      beta[i,] <- (dens[i+1,]*beta[i+1,]) %*% trDens[i,,]
      beta[i,] <- beta[i,]/sum(beta[i,])
    }
  }
  
  for (case in 1 : lt) {
    initP <- initDens * dens[starting[case],] * beta[starting[case],]
    states[starting[case]] <- which(rmultinom(1,1,initP) == 1)
    if (ending[case] == starting[case]) {
      # chain with only 1 states
      next
    }
    for (i in (starting[case]+1) : ending[case]) {
      statesP <- trDens[i-1,states[i-1],] * dens[i,] * beta[i,]
      states[i] <- which(rmultinom(1,1,statesP) == 1)
    }
  }
  
  states
}

NIWSamplerP <- function(data, states, mu0, K0, v0, A0, pd, gap, mu_old=NULL, sigma_old=NULL, lvl=0.8) {
  k <- nrow(mu0)
  mu_all <- matrix(0, k, 2)
  sigma_all <- array(rep(c(1,0,0,1),k), dim = c(2,2,k))
  
  ###############################
  for (kk in c(1,2)) {
    # For empty cluster
    if(sum(states == kk) <= 1) {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
      mu_all[kk,] <- mu_old[kk,]# rnorm(1, mu0[kk,], sigma_all[,,kk]/K0[kk])
      # cat("strange: ", kk, "\n")
      next
    } 
    # First compute the mu of nDML
    # only for 1,2,5
    X <- data[states == kk,]
    X <- reformData(X, pd)
    # 3 indexes, i1: both observed, i2: X only, i3: Y only
    i1 <- which(!is.na(X[,1]) & !is.na(X[,2]))
    i2 <- which(!is.na(X[,1]) & is.na(X[,2]))
    i3 <- which(is.na(X[,1]) & !is.na(X[,2]))
    n <- c(length(i1),length(i2),length(i3))
    
    X_bar <- mean(X[c(i1,i2),1])
    Y_bar <- mean(X[c(i3),2])
    sq0 <- K0[kk] * (sigma_all[1,1,kk]+sigma_all[2,2,kk])
    sq1 <- (n[1]+n[2]) * (sigma_all[1,1,kk]+sigma_all[2,2,kk])
    sq2 <- n[3] * sigma_all[1,1,kk]
    
    mu_n <- (mu0[kk,1]*sq0 + X_bar*sq1 + Y_bar*sq2) / (sq0+sq1+sq2)
    sig_n <- sigma_all[1,1,kk] * (sigma_all[1,1,kk]+sigma_all[2,2,kk]) / (sq0+sq1+sq2)
    
    mu_all[kk,1] <- rnorm(1, mu_n, sqrt(sig_n))
    rm(X)
  }    
  
  ###############################
  kk = 1
  acc <- -1
  oldLL <- sigmaPLL(data, mu_all, sigma_old, v0, A0, pd, states)
  temp <- rnorm(1,sqrt(sigma_old[2,2,kk]), 0.025)
  sigma_all[2,2,kk] <- temp^2
  for (i in c(1,2)) {
    temp <- rnorm(1,sqrt(sigma_old[1,1,i]), 0.025)
    sigma_all[1,1,i] <- temp^2
    sigma_all[2,2,i] <- sigma_all[2,2,1]
  }
  
  newLL <- sigmaPLL(data, mu_all, sigma_all, v0, A0, pd, states)
  
  acc <- exp(newLL - oldLL)
  if (acc < runif(1)) {
    for (i in c(1,2)) {
      sigma_all[,,i] <- sigma_old[,,i]
    }
  }
  rm(acc)
  #############################
  
  for (kk in c(3,4)) {
    # For empty cluster
    if(sum(states == kk) <= 1) {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
      mu_all[kk,] <- mu_old[kk,]# rnorm(1, mu0[kk,], sigma_all[,,kk]/K0[kk])
      # cat("strange: ", kk, "\n")
      next
    } 
    
    # compute the SIGMA of DML
    # only for 3,4
    X <- data[states == kk,]
    X <- reformData(X, pd)
    
    # compute the mu of DML
    if (is.null(mu_old)) {
      i1 <- which(!is.na(X[,1]) & !is.na(X[,2]))
      n <- length(i1)
      mu <- colMeans(X, na.rm = TRUE)
      S <- var(X[i1,]) * (n-1)
      Kn <- K0[kk] + n
      vn <- v0[kk] + n
      mu_n <- (K0[kk] * mu0[kk,] + n * mu) / Kn
      mu_all[kk,] <- rmvnorm(1, mu_n, sigma_old[,,kk]/Kn)
      
    } else {
      # Truncated Normal Prior
      # Use MH sampling      
      
      acc <- -1
      mu_all[kk,] <- rmvnorm(1, mu_old[kk,], diag(2)*0.1)
      
      oldLL <- sum(dm_dmvnorm(X, mu_old[kk,], sigma_old[,,kk], TRUE)) + 
        dmvnorm(mu_old[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
      
      if (kk == 4) {
        #Hypo DMR
        if ((mu_all[kk,1] - mu_all[kk,2]) > gap[1]) {
          newLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) +
            dmvnorm(mu_all[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
          acc <- exp(newLL - oldLL)
        } else {
          if ((mu_old[kk,1] - mu_old[kk,2]) < gap[1]) {
            acc <- 0.5
          } else {
            acc <- -1
          }
        }
      } else if (kk == 3){
        #Hyper DMR
        if ((mu_all[kk,1] - mu_all[kk,2]) < -gap[2]) {
          newLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) +
            dmvnorm(mu_all[kk,], mu0[kk,], diag(2)/K0[kk], TRUE)
          acc <- exp(newLL - oldLL)
        } else {
          if ((mu_old[kk,1] - mu_old[kk,2]) > -gap[2]) {
            acc <- 0.5
          } else {
            acc <- -1
          }
        }
      }
      
      if (acc < runif(1)) {
        mu_all[kk,] <- mu_old[kk,]
      }
    }
    
    oldLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_old[,,kk], TRUE)) + 
      log(diwish(sigma_old[,,kk], v0[kk], A0[,,kk]))
    
    s1_old <- sqrt(sigma_old[1,1,kk])
    s2_old <- sqrt(sigma_old[2,2,kk])
    #rho_old <- tan(sigma_old[1,2,kk]/abs(s1_old*s2_old)*pi/2)
    rho_old <- sigma_old[1,2,kk]/abs(s1_old*s2_old)
    
    s1 <- rnorm(1, s1_old, 0.025)
    s2 <- rnorm(1, s2_old, 0.025)
    rho <- rnorm(1, rho_old, 0.01)
    while (rho <= -1 | rho >= 1) {
      rho <- rnorm(1, rho_old, 0.025)
    }
    sigma_all[1,1,kk] <- s1^2
    sigma_all[2,2,kk] <- s2^2
    #sigma_all[1,2,kk] <- abs(s1*s2)*tanh(rho)*2/pi
    sigma_all[1,2,kk] <- abs(s1*s2)*rho
    sigma_all[2,1,kk] <- sigma_all[1,2,kk]
    
    lvl = 0.8
    pts <- ellipse(sigma_all[,,kk], centre = mu_all[kk,], level = lvl)
    dd <- pts[,2] - pts[,1]
    # New cov is ok
    if (all(dd > 0) | all(dd < 0)) {
      newLL <- sum(dm_dmvnorm(X, mu_all[kk,], sigma_all[,,kk], TRUE)) + 
        log(diwish(sigma_all[,,kk], v0[kk], A0[,,kk]))
      acc <- exp(newLL - oldLL)
      if (acc >= runif(1)) {
        next
      }
    }
    
    # New cov is not ok, check if old cov is ok
    pts <- ellipse(sigma_old[,,kk], centre = mu_all[kk,], level = lvl)
    dd <- pts[,2] - pts[,1]
    if (all(dd > 0) | all(dd < 0)) {
      sigma_all[,,kk] <- sigma_old[,,kk]
    } else {
      sigma_all[,,kk] <- riwish(v0[kk], A0[,,kk])
    }
    rm(X)
  }
  
  return(list(mu = mu_all, sigma = sigma_all))
}

sigmaPLL <- function(data, mu, sigma, v0, A0, pd, states) {
  
  # LL from prior
  LL <- log(dinvgamma(sigma[2,2,1], A0[2,2,2], A0[2,2,1]))
  for (kk in c(1,2)) {
    if(!any(states == kk)) {
      next
    }
    # prior
    LL <- LL + log(dinvgamma(sigma[1,1,kk], v0[kk], A0[1,1,kk]))
    
    X <- data[states == kk,]
    X <- reformData(X, pd)
    i1 <- which(!is.na(X[,1]) & !is.na(X[,2]))
    i2 <- which(!is.na(X[,1]) & is.na(X[,2]))
    i3 <- which(is.na(X[,1]) & !is.na(X[,2]))
    
    LL <- LL + sum(dnorm(X[i1,1], mu[kk,1], sqrt(sigma[1,1,kk]), log = TRUE))
    LL <- LL + sum(dnorm(X[i1,2]-X[i1,1], 0, sqrt(sigma[2,2,kk]), log = TRUE))
    
    LL <- LL + sum(dnorm(X[i2,1], mu[kk,1], sqrt(sigma[1,1,kk]), log = TRUE))
    LL <- LL + sum(dnorm(X[i3,2], mu[kk,1], sqrt(sigma[1,1,kk]+sigma[2,2,kk]), log = TRUE))
  }
  LL
}

Ssampler <- function(logDens, initDens = rep(1/ncol(logDens),ncol(logDens))) {
  const <- apply(logDens,1,max)
  logDens <-  logDens - const
  dens <- exp(logDens)
  
  t_m <- nrow(logDens)
  k <- ncol(logDens)
  
  states <- rep(0, t_m)
  for (i in 1 : t_m) {
    states[i] <- which.max(rmultinom(1,1,dens[i,]*initDens))
  }
  states
}