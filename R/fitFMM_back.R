###############################################################
# Internal function: fit multicomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   nback: number of FMM components to be fitted.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   objectFMM: FMM object to refine the fitting process.
#   staticComponents: fixed components of previous objectFMM.
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
                      objectFMM = NULL, staticComponents = NULL,
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
    cat("|\n")
    cat("|")
    completedPercentage <- 0.00001
    previousPercentage <- completedPercentage
  }

  # Object initialization.
  predichosComponente <- list()
  ajusteComponente <- list()

  for(i in 1:nback){
    predichosComponente[[i]] <- rep(0,n)
  }
  prevAdjMob <- NULL

  # Backfitting algorithm: iteration
  for(i in 1:maxiter){

    # Backfitting algorithm: component
    for(j in 1:nback){

      # data for component j: difference between vData and all other components fitted values
      vDataAjuste <- vData
      for(k in 1:nback){
        if(j != k) vDataAjuste <- vDataAjuste - predichosComponente[[k]]
      }

      # component j fitting using fitFMM_unit function
      ajusteComponente[[j]] <- fitFMM_unit(vDataAjuste,timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                            lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid[[j]], omegaMax = omegaMax,
                                            omegaGrid = omegaGrid[[j]], numReps = numReps, usedApply)
      predichosComponente[[j]] <- ajusteComponente[[j]]@fittedValues
      # showProgress
      if(showProgress){
        completedPercentage <- completedPercentage + 100/(nback*maxiter)
        if(ceiling(previousPercentage) < floor(completedPercentage)){
          numMarcas <- sum((seq(ceiling(previousPercentage), floor(completedPercentage), by=1)%%partialMarkLength == 0))
        } else {
          numMarcas <- 0
        }
        if (numMarcas > 0) {
          for(m in 1:numMarcas) cat("=")
          previousPercentage <- completedPercentage
        }
      }

    }

    # Check stop criterion
    # Fitted values as sum of all components
    adjMob <- rep(0,n)
    for(j in 1:nback){
      adjMob <- adjMob + predichosComponente[[j]]
    }
    if(!is.null(prevAdjMob)){

      if(PV(vData, prevAdjMob) > PV(vData, adjMob)){
        ajusteComponente <- ajusteComponenteAnt
        adjMob <- prevAdjMob
        break
      }

      if(stopFunction(vData,adjMob,prevAdjMob)){
        break
      }

    }

    prevAdjMob <- adjMob
    ajusteComponenteAnt <- ajusteComponente

  }

  # showProgress
  if(showProgress){
    if(completedPercentage < 100){
      completedPercentage <- 100
      if(ceiling(previousPercentage) < floor(completedPercentage)){
        numMarcas <- sum((seq(ceiling(previousPercentage),floor(completedPercentage),by=1)%%partialMarkLength == 0))
      } else {
        numMarcas <- 0
      }
      if (numMarcas > 0) {
        for(m in 1:numMarcas) cat("=")
        previousPercentage <- completedPercentage
      }
    }
    cat("|\n")
    if(i == maxiter){
      cat("Stopped by reaching maximum iterations\n")
    } else {
      cat("Stopped by the stopFunction\n")
    }
  }

  # alpha, beta y omega estimates
  alpha <- rep(0,nback)
  beta <- rep(0,nback)
  omega <- rep(0,nback)
  for(j in 1:nback){
    if(j %in% staticComponents){
      alpha[j] <-getAlpha(objectFMM)[j]
      beta[j] <- getBeta(objectFMM)[j]
      omega[j] <- getOmega(objectFMM)[j]
    } else {
      alpha[j] <- getAlpha(ajusteComponente[[j]])
      beta[j] <- getBeta(ajusteComponente[[j]])
      omega[j] <- getOmega(ajusteComponente[[j]])
    }
  }

  # A and M estimates are recalculated by linear regression
  cos.phi <- list()
  for(j in 1:nback){
    cos.phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints - alpha[j])/2)))
  }
  M <- matrix(unlist(cos.phi),ncol=nback)
  regresion <- lm(vData ~ M)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]

  # Fitted values
  adjMob <- predict(regresion)

  # Residual sum of squares
  SSE <- sum((adjMob - vData)^2)

  names(A) <- paste("A", 1:length(A), sep="")

  # Returns an object of class FMM.
  outMobius <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues= adjMob,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega)
  )

  return(outMobius)

}
