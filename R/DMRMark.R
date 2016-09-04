DMRMark <- function(mv, L=rep(1,nrow(mv)), starting=NULL, pd=NULL,
                    initHeuristic=TRUE,
                    GSoptions=NULL) {
  
  # initHeuristic:Heuristic for faster computation.
  #               Much faster. Using heuristic for finding good initial value
  #               and quick convergence. Mask GS controls of GSoptions.
  #               Generally performs good, but NO guarantee of convergence.
  #               Recommanded for starting new anaylsis, getting quick insight.
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
  
  opts <- NULL
  if (is.null(GSoptions)) {
    mv2 <- reformData(mv, pd)
    opts <- MakeGSoptions(D_mu = boundFinder(mv2, 0.1))
    rm(mv2)
  } else {
    opts <- GSoptions
  }
  opts <- oldOption(opts)

  print("Start...")
  if (initHeuristic) {    
    pars <- myHeuristic(mv, pd, L, starting, opts)
    pars <- toNewPars(pars)
    print("Parameter Estimation finished.")
    return(pars)
  }
  
  pars <- FullSample(mv, L, starting, opts$alpha, opts$mu0, opts$K0, opts$v0, opts$A0, pd, opts$cmu0, 
                       gap=opts$gap, lvl=opts$lvl, 
                       burnin = opts$burnin, sampling = opts$sampling, 
                       track = opts$track, withHMM = opts$onHMM, 
                       verbose = opts$verbose, nsamples = opts$nsamples, sampleSep = opts$sampleSep)
  
  if (!opts$track) {
    pars <- toNewPars(pars)
  }
  print("Parameter Estimation finished.")
  return(pars)
}
                          