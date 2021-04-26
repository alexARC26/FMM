################################################################################
# Auxiliary internal functions
# Functions:
#   step1FMM:     M, A and beta initial parameter estimations.
#   bestStep1:    to find the optimal initial parameters estimation.
#   step2FMM:     second step of FMM fitting process.
#   refineFMM:    fitFMM from a previous objectFMM.
#   PV:           percentage of variability explained.
#   PVj:          percentage of variability explained by each component of
#                 FMM model.
#   angularmean:  to compute the angular mean.
#   seqTimes:     to build a sequence of equally time points spaced in range
#                 [0,2*pi].
################################################################################


################################################################################
# Internal function: to estimate M, A and beta initial parameters
# also returns residual sum of squared (RSS).
# Arguments:
#    alphaOmegaParameters: vector of the parameters alpha and omega
#    vData: data to be fitted an FMM model.
#    timePoints: one single period time points.
# Returns a 6-length numerical vector: M, A, alpha, beta, omega and RSS
################################################################################
step1FMM <- function(alphaOmegaParameters, vData, timePoints) {

  alphaParameter <- alphaOmegaParameters[1]
  omegaParameter <- alphaOmegaParameters[2]

  mobiusTerm <- 2*atan(omegaParameter*tan((timePoints - alphaParameter)/2))
  tStar <- alphaParameter + mobiusTerm

  # Given alpha and omega, a cosinor model is computed with t* in
  # order to obtain delta (cosCoeff) and gamma (sinCoeff).
  # Linear Model exact expressions are used to improve performance.
  costStar <- cos(tStar)
  sentstar <- sin(tStar)
  covMatrix <- stats::cov(cbind(vData, costStar, sentstar))
  denominator <- covMatrix[2,2]*covMatrix[3,3] - covMatrix[2,3]^2
  cosCoeff <- (covMatrix[1,2]*covMatrix[3,3] -
                 covMatrix[1,3]*covMatrix[2,3])/denominator
  sinCoeff <- (covMatrix[1,3]*covMatrix[2,2] -
                 covMatrix[1,2]*covMatrix[2,3])/denominator
  mParameter <- mean(vData) - cosCoeff*mean(costStar) - sinCoeff*mean(sentstar)

  phiEst <- atan2(-sinCoeff, cosCoeff)
  aParameter <- sqrt(cosCoeff^2 + sinCoeff^2)
  betaParameter <- (phiEst+alphaParameter)%%(2*pi)

  mobiusRegression <- mParameter + aParameter*cos(betaParameter + mobiusTerm)
  residualSS <- sum((vData - mobiusRegression)^2)/length(timePoints)

  return(c(mParameter, aParameter, alphaParameter, betaParameter,
           omegaParameter, residualSS))
}


################################################################################
# Internal function: to find the optimal initial parameter estimation
# Arguments:
#    vData: data to be fitted an FMM model.
#    step1: a data.frame with estimates of
#           M, A, alpha, beta, omega, RSS as columns.
# Returns the optimal row of step1 argument.
# optimum: minimum RSS with several stability conditions.
################################################################################
bestStep1 <- function(vData, step1){

  # step1 in decreasing order by RSS
  orderedModelParameters <- order(step1[,"RSS"])

  maxVData <- max(vData)
  minVData <- min(vData)
  n <- length(vData)

  # iterative search: go through rows ordered step 1
  #    until the first one that verifies the stability conditions
  bestModelFound <- FALSE
  i <- 1
  while(!bestModelFound){
    # parameters
    mParameter <- step1[orderedModelParameters[i], "M"]
    aParameter <- step1[orderedModelParameters[i], "A"]
    alphaParameter <- step1[orderedModelParameters[i], "alpha"]
    betaParameter <- step1[orderedModelParameters[i], "beta"]
    omegaParameter <- step1[orderedModelParameters[i], "omega"]
    sigma <- sqrt(step1[orderedModelParameters[i], "RSS"]*n/(n-5))

    # stability conditions
    amplitudeUpperBound <- mParameter + aParameter
    amplitudeLowerBound <- mParameter - aParameter
    rest1 <- amplitudeUpperBound <= maxVData + 1.96*sigma
    rest2 <- amplitudeLowerBound >= minVData - 1.96*sigma

    # it is necessary to check that there are no NA,
    # because it can be an extreme solution
    if(is.na(rest1)) rest1 <- FALSE
    if(is.na(rest2)) rest2 <- FALSE

    if(rest1 & rest2){
      bestModelFound <- TRUE
    } else {
      i <- i+1
    }

    if(i > nrow(step1))
      return(NULL)
  }

  return(step1[orderedModelParameters[i],])
}

################################################################################
# Internal function: second step of FMM fitting process
# Arguments:
#   parameters: M, A, alpha, beta, omega initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omegaMax: max value for omega.
################################################################################
step2FMM <- function(parameters, vData, timePoints, omegaMax){

  n <- length(timePoints)

  # FMM model and residual sum of squares
  modelFMM <- parameters[1] + parameters[2] *
    cos(parameters[4]+2*atan2(parameters[5]*sin((timePoints - parameters[3])/2),
                                cos((timePoints - parameters[3])/2)))
  residualSS <- sum((modelFMM - vData)^2)/n
  sigma <- sqrt(residualSS*n/(n - 5))

  # When amplitude condition is valid, it returns RSS
  # else it returns infinite.
  amplitudeUpperBound <- parameters[1] + parameters[2]
  amplitudeLowerBound <- parameters[1] - parameters[2]
  rest1 <- amplitudeUpperBound <= max(vData) + 1.96*sigma
  rest2 <- amplitudeLowerBound >= min(vData) - 1.96*sigma

  # Other integrity conditions that must be met
  rest3 <- parameters[2] > 0  # A > 0
  rest4 <- parameters[5] > 0  &  parameters[5] <= omegaMax # omega > 0 and omega <= omegaMax
  if(rest1 & rest2 & rest3 & rest4)
    return(residualSS)
  else
    return(Inf)
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   the FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values.
################################################################################
PV <- function(vData,pred){
  meanVData <- mean(vData)
  return(1 - sum((vData-pred)^2)/sum((vData-meanVData)^2))
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   each component of FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   alpha, beta, omega: vectors of corresponding parameter estimates.
################################################################################
PVj <- function(vData, timePoints, alpha, beta, omega){

  # fitted values of each wave
  nComponents <- length(alpha)
  w <- list()
  for(i in 1:nComponents){
    w[[i]] <- cos(beta[i] + 2*atan(omega[i]*tan((timePoints-alpha[i])/2)))
  }

  # The fitting is recalculated only up to wave i and
  # the percentage of variability explained is determined
  PV_hasta <- c()
  for(i in 1:nComponents){
    M <- matrix(unlist(w[1:i]),ncol=i)
    regresion <- lm(vData ~ M)
    predichos <- predict(regresion)
    PV_hasta[i] <- PV(vData,predichos)
  }

  # individual percentage of variability is the part that adds to the whole
  PV_individual <- PV_hasta
  for(i in 2:nComponents){
    PV_individual[i] <- PV_hasta[i]-sum(PV_individual[1:(i-1)])
  }

  return(PV_individual)

}

################################################################################
# Internal function: to build a sequence of equally time points spaced
#                    in range [0,2*pi].
# Arguments:
#   n: secuence length.
################################################################################
seqTimes <- function(n){
  timePoints<-seq(0,2*pi,by=2*pi/n)
  timePoints<-timePoints[-length(timePoints)]
  return(timePoints)
}

################################################################################
# Internal function: to compute the angular mean.
# Arguments:
#   angles: input vector of angles.
################################################################################
#Computes the Angular Mean
angularmean <- function(angles){
  n <- length(angles)
  a.mean <- atan2(sum(sin(angles)),sum(cos(angles)))
  return(a.mean)
}

################################################################################
# Internal function: to replicate a grid and return a list with replications.
# Arguments:
#   angles: grid to replicate.
#   nback: times the grid is going to be replicated
################################################################################
replicateGrid <- function(grid, nback){
  return(replicate(n = nback, grid, simplify = FALSE))
}

################################################################################
# Internal function: to calculate cosine term of FMM model.
# Arguments:
#   alpha, beta, omega: parameters
#   timePoints: time poinst where FMM model is computed
################################################################################
calculateCosPhi <- function(alpha, beta, omega, timePoints){
  calculateSingleCosPhi <- function(alpha, beta, omega){
    return(cos(beta + 2*atan(omega*tan((timePoints - alpha)/2))))
  }
  return(mapply(FUN = calculateSingleCosPhi, alpha = alpha, beta = beta, omega = omega))
}

################################################################################
# Internal function: return parallelized apply function depending on the OS.
# Returns function to be used.
################################################################################
getApply <- function(parallelize = FALSE, nCores = min(12, parallel::detectCores() - 1)){

  getApply_Rbase <- function(){
    usedApply <- function(FUN, X, ...) t(apply(X = X, MARGIN = 1, FUN = FUN, ...))
  }

  getParallelApply_Windows <- function(parallelCluster){
    usedApply <- function(FUN, X, ...) t(parallel::parApply(parallelCluster, FUN = FUN,
                                                            X = X, MARGIN = 1, ...))
    return(usedApply)
  }

  parallelFunction_Unix<-function(nCores){
    # A parallelized apply function does not exist, so it must be translated to a lapply
    usedApply <- function(FUN, X, ...){
      matrix(unlist(parallel::mclapply(X = asplit(X, 1), FUN = FUN, mc.cores = nCores, ...)),
             nrow = nrow(X), byrow = T)
    }
    return(usedApply)
  }

  if(parallelize){
    # different ways to implement parallelization depending on OS:
    if(.Platform$OS.type == "windows"){
      parallelCluster <- parallel::makePSOCKcluster(nCores)
      doParallel::registerDoParallel(parallelCluster)
      usedApply <- getParallelApply_Windows(parallelCluster)
    }else{
      usedApply <- parallelFunction_Unix(nCores)
      parallelCluster <- NULL
    }
  }else{
    # R base apply:
    usedApply <- getApply_Rbase()
    parallelCluster <- NULL
  }

  return(list(usedApply, parallelCluster))
}

################################################################################
#                      RESTRICTED FMM INTERNAL FUNCTIONS                       #
################################################################################

################################################################################
# Internal function: iterate fitting on omegasIter grid (FMM_unit_restr models
# with fixed omega)
# Arguments:
#   vData: data to be fitted an FMM model.
#   omegasIter: omegas grid (depending on original omega grid and on
#               restrictions)
#   betaRestrictions: beta's constraint vector.
#   timePoints: one single period time points.
#   lengthAlphaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid: grid of alpha parameter.
#   numReps: number of times the alpha-omega grid search is repeated.
#   nback: number of FMM components to fit.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   parallelize: TRUE to use parallelized procedure to fit FMM model
# Returns a list with objects of class FMM.
################################################################################
iterateOmegaGrid <- function(vData, omegasIter, betaRestrictions, timePoints = seqTimes(length(vData)),
                       lengthAlphaGrid = 48, alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid),
                       numReps = numReps, nback = nback, maxiter = maxiter, parallelize = FALSE){
  if(parallelize){
    # cores for parallelization registered: 'foreach' is used
    objectFMMList <- foreach::foreach(omegas = iterators::iter(omegasIter, by="row")) %dopar% {
      backFittingRestr(vData = vData, timePoints = timePoints, omegas = omegas,
                       lengthAlphaGrid = lengthAlphaGrid, alphaGrid = alphaGrid, numReps = numReps,
                       nback = nback, maxiter = maxiter, betaRestrictions = betaRestrictions)
    }
  }else{
    # 'for' loop (non-parallelized version)
    objectFMMList <- list()
    for(i in 1:length(omegasIter[,1])){
      omegas <- omegasIter[i,]
      suppressWarnings(
        objectFMMList[[i]] <-  backFittingRestr(vData = vData, timePoints = timePoints,
                                                omegas = omegas, lengthAlphaGrid = lengthAlphaGrid,
                                                alphaGrid = alphaGrid, numReps = numReps, nback = nback,
                                                maxiter = maxiter, betaRestrictions = betaRestrictions)
      )
    }
  }
  return(objectFMMList)

}


#######################################################################################
# Internal function: backfitting algorithm with fixed omegas.
# Arguments:
#   vData: data to be fitted an FMM model.
#   omegas: fixed omegas to fit FMM models.
#   betaRestrictions: beta's constraint vector.
#   timePoints: one single period time points.
#   alphaGrid: grid of alpha parameter.
#   numReps: number of times the alpha-omega grid search is repeated.
#   nback: number of FMM components to fit.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   parallelize: TRUE to use parallelized procedure to fit FMM model
#   stopFunction: Function to check the convergence criterion for the
#                 backfitting algorithm
# Returns an object of class FMM.
#######################################################################################
backFittingRestr <- function(vData, omegas, nback, betaRestrictions,
                             timePoints = seqTimes(length(vData)), lengthAlphaGrid = 48,
                             alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid),
                             numReps = 3, maxiter = nback, stopFunction = alwaysFalse){
  n <- length(vData)
  omegas <- as.numeric(omegas)
  fittedValuesPerComponent <- matrix(0, ncol = nback, nrow = n)
  fittedFMMPerComponent <- list()
  prevFittedFMMvalues <- NULL

  # Backfitting algorithm: iteration
  for(i in 1:maxiter){
    # Backfitting algorithm: component
    for(j in 1:nback){
      # data for component j: difference between vData and all other components fitted values
      backFittingData <- vData - apply(as.matrix(fittedValuesPerComponent[,-j]), 1, sum)
      # component j fitting using fitFMM_unit_restr function
      fittedFMMPerComponent[[j]] <- fitFMM_unit_restr(backFittingData, omegas[j], timePoints = timePoints,
                                                      lengthAlphaGrid = lengthAlphaGrid, alphaGrid = alphaGrid[[j]],
                                                      numReps = numReps)
      fittedValuesPerComponent[,j] <- fittedFMMPerComponent[[j]]@fittedValues
    }

    # Check stop criterion
    # Fitted values as sum of all components
    fittedFMMvalues <- apply(fittedValuesPerComponent, 1, sum)

    if(!is.null(prevFittedFMMvalues)){

      if(PV(vData, prevFittedFMMvalues) > PV(vData, fittedFMMvalues)){
        fittedFMMPerComponent <- previousFittedFMMPerComponent
        fittedFMMvalues <- prevFittedFMMvalues
        break
      }
      if(stopFunction(vData, fittedFMMvalues, prevFittedFMMvalues)){
        break
      }
    }
    prevFittedFMMvalues <- fittedFMMvalues
    previousFittedFMMPerComponent <- fittedFMMPerComponent
  }

  # alpha, beta y omega estimates
  alpha <- unlist(lapply(fittedFMMPerComponent, getAlpha))
  beta <- unlist(lapply(fittedFMMPerComponent, getBeta))
  omega <- unlist(lapply(fittedFMMPerComponent, getOmega))

  # beta restrictions: calculate angular mean of beta parameters
  # the nearest betas are chosen
  restBeta <- beta
  markedBetas <- rep(0, length(betaRestrictions))
  betaIndexVector <- 1:length(betaRestrictions)
  for(indRes in unique(betaRestrictions)){
    numComponents <- sum(betaRestrictions == indRes)
    betaIndex <- betaIndexVector[markedBetas == 0][1]
    distance <- abs(restBeta - restBeta[betaIndex])
    distance[markedBetas == 1] <- Inf
    nearestBetasIndex <- order(distance)[1:numComponents]
    markedBetas[nearestBetasIndex] <- 1
    restBeta[nearestBetasIndex] <- angularmean(restBeta[nearestBetasIndex])%%(2*pi)
  }

  # A and M estimates are recalculated by linear regression
  cosPhi <- calculateCosPhi(alpha = alpha, beta = restBeta, omega = omega, timePoints = timePoints)
  regresion <- lm(vData ~ cosPhi)

  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]
  fittedFMMvalues <- predict(regresion)

  # Residual sum of squares
  SSE <- ifelse(sum(A < 0, na.rm = TRUE) > 0, Inf, sum((fittedFMMvalues-vData)^2))
  names(A) <- paste("A", 1:length(A), sep="")

  # Returns an object of class FMM
  outMobius <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = restBeta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = i
  )

  return(outMobius)
}

#######################################################################################
# Internal function: to fit monocomponent FMM models with fixed omega.
# Arguments:
#   vData: data to be fitted an FMM model.
#   omega: value of the omega parameter.
#   timePoints: one single period time points.
#   lengthAlphaGrid: precision of the grid of alpha parameter.
#   alphaGrid: grid of alpha parameter.
#   numReps: number of times the alpha-omega grid search is repeated.
# Returns an object of class FMM.
#######################################################################################
fitFMM_unit_restr<-function(vData, omega, timePoints = seqTimes(length(vData)),
                            lengthAlphaGrid = 48, alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                            numReps = 3, usedApply = getApply(FALSE)[[1]]){
  n <- length(vData)

  ## Step 1: initial values of M, A, alpha, beta and omega
  # alpha and omega are fixed and cosinor model is used to calculate the
  # rest of the parameters.
  # step1FMM function is used to make this estimate
  grid <- expand.grid(alphaGrid, omega)
  step1 <- usedApply(FUN = step1FMM, X = grid, vData = vData, timePoints = timePoints)
  colnames(step1) <- c("M","A","alpha","beta","omega","RSS")

  # We find the optimal initial parameters,
  # minimizing Residual Sum of Squared with several stability conditions.
  # We use bestStep1 internal function
  bestPar <- bestStep1(vData,step1)

  # When the fixed omega is so extreme that the fitting is no possible,
  # the null fitted model is returned.
  if(is.null(bestPar)){
    outMobius <- FMM(
      M = 0,
      A = 0,
      alpha = 0,
      beta = 0,
      omega = omega,
      timePoints = timePoints,
      summarizedData = vData,
      fittedValues = rep(0,length(vData)),
      SSE = sum(vData^2),
      R2 = PV(vData, rep(0,length(vData))),
      nIter = 0
    )
    return(outMobius)
  }

  ## Step 2: Nelder-Mead optimization. 'step2FMM_restr' function is used.
  if(!is.infinite(step2FMM_restr(bestPar[1:4], vData = vData, timePoints = timePoints, omega = omega))){
    nelderMead <- optim(par = bestPar[1:4], fn = step2FMM_restr, vData = vData,
                        timePoints = timePoints, omega = omega)
    parFinal <- c(nelderMead$par, omega)
  } else {
    parFinal <- c(bestPar[1:4],omega)
  }

  names(parFinal) <- c("M", "A", "alpha", "beta", "omega")
  fittedFMMvalues <- parFinal["M"] + parFinal["A"]*
    cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((fittedFMMvalues-vData)^2)

  # alpha and beta between 0 and 2pi
  parFinal[3] <- parFinal[3]%%(2*pi)
  parFinal[4] <- parFinal[4]%%(2*pi)

  # the grid search is repeated numReps
  numReps <- numReps - 1
  while(numReps > 0){
    # new grid for alpha between 0 and 2pi
    nAlphaGrid <- length(alphaGrid)
    amplitudeAlphaGrid <- 1.5*mean(diff(alphaGrid))
    alphaGrid <- seq(parFinal[3]-amplitudeAlphaGrid,parFinal[3]+amplitudeAlphaGrid,
                     length.out = nAlphaGrid)
    alphaGrid <- alphaGrid%%(2*pi)

    ## Step 1: initial parameters
    grid <- as.matrix(expand.grid(alphaGrid,omega))
    step1 <- t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints))
    colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
    antBestPar <- bestPar
    bestPar <- bestStep1(vData,step1)

    # None satisfies the conditions
    if(is.null(bestPar)){
      bestPar <- antBestPar
      numReps <- 0
      warning("FMM model may be no appropiate")
    }

    ## Step 2: Nelder-Mead optimization
    if(!is.infinite(step2FMM_restr(bestPar[1:4],vData = vData, timePoints = timePoints, omega = omega))){
      nelderMead <- optim(par = bestPar[1:4], fn = step2FMM_restr, vData = vData,
                          timePoints = timePoints, omega = omega)
      parFinal <- c(nelderMead$par,omega)
    } else {
      parFinal <- c(bestPar[1:4],omega)
    }

    # alpha and beta between 0 and 2pi
    parFinal[3] <- parFinal[3]%%(2*pi)
    parFinal[4] <- parFinal[4]%%(2*pi)

    numReps <- numReps - 1
  }
  names(parFinal) <- c("M","A","alpha","beta","omega")

  # Returns an object of class FMM
  fittedFMMvalues <- parFinal["M"] + parFinal["A"]*
    cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((fittedFMMvalues-vData)^2)

  outMobius <- FMM(
    M = parFinal["M"],
    A = parFinal["A"],
    alpha = parFinal[3],
    beta = parFinal[4],
    omega = parFinal[5],
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = 0,
    nIter = 0
  )
  return(outMobius)
}

################################################################################
# Internal function: to optimize omega.
# It is used in the extra optimization step of omega,
# within fitFMM_restr function.
# Arguments:
#   uniqueOmegas: grid of omega parameters.
#   indOmegas: omega's constraint vector.
#   objFMM: FMM object to refine the omega fitting.
#   omegaMax: max value for omega.
# Returns the residual sum of squares or
# 'Inf' when the integrity conditions are not met.
################################################################################
stepOmega <- function(uniqueOmegas, indOmegas, objFMM, omegaMax){

  timePoints <- getTimePoints(objFMM)
  n <- length(timePoints)
  M <- getM(objFMM)
  A <- getA(objFMM)
  alpha <- getAlpha(objFMM)
  beta <- getBeta(objFMM)
  omega <- uniqueOmegas[indOmegas]
  nback <- length(alpha)
  vData <- getSummarizedData(objFMM)

  # FMM fitting and RSS
  fittedValuesPerComponent <- matrix(0, ncol = nback, nrow = n)
  for(j in 1:nback){
    fittedValuesPerComponent[,j] <- A[j]*cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  fittedValues <- apply(fittedValuesPerComponent, 1, sum) + M
  RSS <- sum((fittedValues - vData)^2)

  # Other integrity conditions that must be met
  rest1 <- all(uniqueOmegas <= omegaMax)
  rest2 <- all(uniqueOmegas >= 0)
  if(rest1 & rest2){
    return(RSS)
  }else{
    return(Inf)
  }
}

################################################################################
# Internal function: second step of FMM fitting process with fixed omega
# Arguments:
#   param: M, A, alpha, beta initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omega: fixed value of omega.
################################################################################
step2FMM_restr <- function(parameters, vData, timePoints, omega){

  n <- length(timePoints)
  # FMM model and residual sum of squares
  modelFMM <- parameters[1] + parameters[2] *
    cos(parameters[4]+2*atan2(omega*sin((timePoints - parameters[3])/2),
                              cos((timePoints - parameters[3])/2)))
  residualSS <- sum((modelFMM - vData)^2)/n
  sigma <- sqrt(residualSS*n/(n - 5))

  # When amplitude condition is valid, it returns RSS
  # else it returns infinite.
  amplitudeUpperBound <- parameters[1] + parameters[2]
  amplitudeLowerBound <- parameters[1] - parameters[2]
  rest1 <- amplitudeUpperBound <= max(vData) + 1.96*sigma
  rest2 <- amplitudeLowerBound >= min(vData) - 1.96*sigma

  # Other integrity conditions that must be met
  rest3 <- parameters[2] > 0 #A > 0
  if(rest1 & rest2 & rest3){
    return(residualSS)
  }else{
    return(Inf)
  }
}
