################################################################################
# Internal function: fit monocomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of
#                                     alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
# Returns an object of class FMM.
################################################################################
fitFMM_unit<-function(vData, timePoints = seqTimes(length(vData)),
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMax = 1,
                      omegaGrid = exp(seq(log(0.0001), log(omegaMax),
                                          length.out = lengthOmegaGrid)),
                      numReps = 3, parallelize = FALSE, useRcpp = FALSE){

  n <- length(vData)
  grid <- expand.grid(alphaGrid,omegaGrid)

  ## Step 1: initial values of M, A, alpha, beta and omega alpha and omega are
  # fixed and cosinor model is used to calculate the rest of the parameters.
  # step1FMM function is used to make this estimate
  # parallelization o rcpp function can be used

  if(parallelize){
    requireNamespace("doParallel", quietly = TRUE)
    nCores <- detectCores() - 1
    registerDoParallel(cores = nCores)
    parallelCluster <- makeCluster(nCores)
    step1 <- t(parApply(parallelCluster, X = grid, 1, FUN = step1FMM,
                        vData = vData, timePoints = timePoints))
    stopCluster(parallelCluster)
  }else if(useRcpp){
    step1 <- t(apply(grid, 1, FUN = step1FMMrcpp, vData = vData,
                     timePoints = timePoints))
  }else{
    step1 <- t(apply(grid, 1, FUN = step1FMM, vData = vData,
                     timePoints = timePoints))
  }
  colnames(step1) <- c("M","A","alpha","beta","omega","RSS")

  # We find the optimal initial parameters,
  # minimizing Residual Sum of Squared with several stability conditions.
  # We use bestStep1 internal function
  bestPar <- bestStep1(vData, step1)

  ## Step 2: Nelder-Mead optimization. 'step2FMM' function is used.
  nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData,
                      timePoints = timePoints, omegaMax = omegaMax,
                      control=list(warn.1d.NelderMead = FALSE))
  parFinal <- nelderMead$par
  SSE <- nelderMead$value*n

  # alpha and beta between 0 and 2pi
  parFinal[3] <- parFinal[3]%%(2*pi)
  parFinal[4] <- parFinal[4]%%(2*pi)

  # the grid search is repeated numReps
  numReps <- numReps - 1
  while(numReps > 0){

    # new grid for alpha between 0 and 2pi
    nAlphaGrid <- length(alphaGrid)
    amplitudeAlphaGrid <- 1.5*mean(diff(alphaGrid))
    alphaGrid <- seq(parFinal[3] - amplitudeAlphaGrid,
                     parFinal[3] + amplitudeAlphaGrid, length.out = nAlphaGrid)
    alphaGrid <- alphaGrid%%(2*pi)

    # new grid for omega between 0 and omegaMax
    nOmegaGrid <- length(omegaGrid)
    amplitudeOmegaGrid <- 1.5*mean(diff(omegaGrid))
    omegaGrid <- seq(max(parFinal[5] - amplitudeOmegaGrid, 0),
                     min(omegaMax,parFinal[5] + amplitudeOmegaGrid),
                     length.out = nOmegaGrid)
    grid <- as.matrix(expand.grid(alphaGrid,omegaGrid))

    # Step 1: initial parameters
    if(parallelize){
      requireNamespace("doParallel", quietly = TRUE)
      nCores <- detectCores() - 1
      registerDoParallel(cores = nCores)
      parallelCluster <- makeCluster(nCores)
      step1 <- t(parApply(parallelCluster, X = grid, 1, FUN = step1FMM,
                          vData = vData, timePoints = timePoints))
      stopCluster(parallelCluster)
    }else if(useRcpp){
      step1 <- t(apply(grid, 1, FUN = step1FMMrcpp, vData = vData,
                       timePoints = timePoints))
    }else{
      step1 <- t(apply(grid, 1, FUN = step1FMM, vData = vData,
                       timePoints = timePoints))
    }
    colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
    prevBestPar <- bestPar
    bestPar <- bestStep1(vData,step1)

    # None satisfies the conditions
    if(is.null(bestPar)){
      bestPar <- prevBestPar
      numReps <- 0
      warning("FMM model may be no appropiate")
    }

    ## Step 2: Nelder-Mead optimization
    nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData,
                        timePoints = timePoints, omegaMax = omegaMax,
                        control = list(warn.1d.NelderMead = FALSE))
    parFinal <- nelderMead$par

    # alpha and beta between 0 and 2pi
    parFinal[3] <- parFinal[3]%%(2*pi)
    parFinal[4] <- parFinal[4]%%(2*pi)

    numReps <- numReps - 1
  }

  names(parFinal) <- c("M","A","alpha","beta","omega")

  # Returns an object of class FMM.
  adjMob <- parFinal["M"] + parFinal["A"]*
    cos(parFinal["beta"] +
        2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((adjMob-vData)^2)

  outMobius <- FMM(
    M = parFinal["M"],
    A = parFinal["A"],
    alpha = parFinal[3],
    beta = parFinal[4],
    omega = parFinal[5],
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = adjMob,
    SSE = SSE,
    R2 = PV(vData, adjMob)
  )
  return(outMobius)
}

