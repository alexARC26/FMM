###############################################################
# Internal function: fit multicomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   nback: number of FMM components to be fitted.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#                         They can be a list with nback elements, each one for an iteration.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
#   showProgress: TRUE to display a progress indicator on the console.
# Returns an object of class FMM.
###############################################################
fitFMM_back<-function(vData, timePoints = seqTimes(length(vData)), nback,
                      maxiter = nback, stopFunction = alwaysFalse,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMax = 1,
                      omegaGrid = exp(seq(log(0.0001),log(omegaMax),
                                          length.out=lengthOmegaGrid)),
                      numReps = 3, showProgress = TRUE, usedApply){

  n <- length(vData)

  if(!is.list(alphaGrid))
    alphaGrid <- replicateGridAsList(grid = alphaGrid, nback = nback)

  if(!is.list(omegaGrid))
    omegaGrid <- replicateGridAsList(grid = omegaGrid, nback = nback)

  if(showProgress){
    totalMarks <- 50
    partialMarkLength <- 2
    cat("|")
    for(m in 1:totalMarks) cat("-")
    cat("|\n|")
    completedPercentage <- 0.00001
    previousPercentage <- completedPercentage
  }

  # Object initialization.
  fittedValuesPerComponent <- matrix(rep(0, n*nback), ncol = nback)
  fittedFMMPerComponent <- list()

  prevFittedFMMvalues <- NULL

  # Backfitting algorithm: iteration
  for(i in 1:maxiter){
    # Backfitting algorithm: component
    for(j in 1:nback){
      # data for component j: difference between vData and all other components fitted values
      backFittingData <- vData - apply(fittedValuesPerComponent[,-j], 1, sum)

      # component j fitting using fitFMM_unit function
      fittedFMMPerComponent[[j]] <- fitFMM_unit(backFittingData, timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                            lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid[[j]], omegaMax = omegaMax,
                                            omegaGrid = omegaGrid[[j]], numReps = numReps, usedApply)
      fittedValuesPerComponent[,j] <- fittedFMMPerComponent[[j]]@fittedValues
      # showProgress
      if(showProgress){
        completedPercentage <- completedPercentage + 100/(nback*maxiter)
        if(ceiling(previousPercentage) < floor(completedPercentage)){
          numMarcas <- sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1) %% partialMarkLength == 0))
          for(m in 1:numMarcas) cat("=")
          previousPercentage <- completedPercentage
        }
      }
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

  # showProgress
  if(showProgress){
    if(completedPercentage < 100){
      completedPercentage <- 100
      if(ceiling(previousPercentage) < floor(completedPercentage)){
        numMarcas <- sum((seq(ceiling(previousPercentage),
                              floor(completedPercentage), by = 1) %% partialMarkLength == 0))
      } else {
        numMarcas <- 0
      }
      if (numMarcas > 0) {
        for(m in 1:numMarcas) cat("=")
        previousPercentage <- completedPercentage
      }
    }
    if(i == maxiter){
      cat("|\nStopped by reaching maximum iterations\n")
    } else {
      cat("|\nStopped by the stopFunction\n")
    }
  }

  alpha <- unlist(lapply(fittedFMMPerComponent, getAlpha))
  beta <- unlist(lapply(fittedFMMPerComponent, getBeta))
  omega <- unlist(lapply(fittedFMMPerComponent, getOmega))

  # A and M estimates are recalculated by linear regression
  cos.phi <- list()
  for(j in 1:nback)
    cos.phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints - alpha[j])/2)))

  designMatrix <- matrix(unlist(cos.phi), ncol = nback)
  linearModel <- lm(vData ~ designMatrix)
  M <- linearModel$coefficients[1]
  A <- linearModel$coefficients[-1]
  fittedFMMvalues <- predict(linearModel)

  # Residual sum of squares
  SSE <- sum((fittedFMMvalues - vData)^2)

  names(A) <- paste("A", 1:length(A), sep = "")

  # Returns an object of class FMM.
  outMobius <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues= fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega)
  )
  return(outMobius)
}
