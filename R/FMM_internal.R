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

  phiEst <- atan2(-sinCoeff, cosCoeff)   # acrophase (phi)
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

replicateGrid <- function(grid, nback){
  return(replicate(n = nback, grid, simplify = FALSE))
}

calculateCosPhi <- function(alpha, beta, omega, timePoints){
  calculateSingleCosPhi <- function(alpha, beta, omega){
    return(cos(beta + 2*atan(omega*tan((timePoints - alpha)/2))))
  }
  return(mapply(FUN=calculateSingleCosPhi, alpha=alpha, beta=beta, omega=omega))
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

