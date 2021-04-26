###############################################################
# Internal functions to fit restricted FMM models
# Functions:
#   fitFMM_restr:           to fit restricted multicomponent FMM models.
#                           An external grid for omega parameter is used.
#   stepOmega:              to optimize omega.
#   fitFMM_unit_restr:      to fit monocomponent FMM models with fixed omega.
#   step2FMM_restr:         second step of FMM fitting process with fixed omega.
#   fitFMM_restr_omegaBeta: to fit restricted multicomponent FMM models.
#                           Nested backfitting algorithm is used.
###############################################################

###############################################################
# Internal function: to fit restricted multicomponent FMM models.
# An external grid for omega parameter is used for fitting process.
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   nback: number of FMM components to be fitted.
#   betaRestrictions: beta's constraint vector.
#   omegaRestrictions: omega's constraint vector.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   objectFMM: FMM object to refine the fitting process.
#   staticComponents: fixed components of previous objectFMM.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#                         They can be a list with nback elements, each one for an iteration.
#   omegaMin: min value for omega.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
#   parallelize: TRUE to use parallelized procedure to fit restricted FMM model.
# Returns an object of class FMM.
###############################################################
fitFMM_restr<-function(vData, timePoints = seqTimes(length(vData)), nback,
                      betaRestrictions, omegaRestrictions, maxiter = nback,
                      stopFunction = alwaysFalse, objectFMM = NULL, staticComponents = NULL,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMin = 0.0001, omegaMax = 1,
                      omegaGrid = exp(seq(log(omegaMin),log(omegaMax), length.out = lengthOmegaGrid)),
                      numReps = 3, parallelize = FALSE){

  n <- length(vData)
  betaRestrictions <- sort(betaRestrictions)
  omegaRestrictions <- sort(omegaRestrictions)
  alphaGrid <- replicateGrid(alphaGrid, nback = nback)

  # External grid of omega parameters
  numOmegas <- length(unique(omegaRestrictions))
  listOmegas <- replicateGrid(omegaGrid, numOmegas)
  gridOmegas <- expand.grid(listOmegas)
  omegasIter <- gridOmegas[,omegaRestrictions]

  # External loop on the omega grid, setting its value
  objectFMMList <- iterateOmegaGrid(vData = vData, omegasIter = omegasIter, betaRestrictions = betaRestrictions,
                                    timePoints = timePoints, alphaGrid = alphaGrid, numReps = numReps,
                                    nback = nback, maxiter = maxiter, parallelize = parallelize)

  # We keep the best solution
  SSElist <- lapply(objectFMMList, getSSE)
  outMobius <- objectFMMList[[which.min(SSElist)]]

  # Extra optimization to allow omega's to move more freely.
  # The stepOmega function is used.
  uniqueOmegas <- unique(getOmega(outMobius))
  uniqueOmegasOptim <- optim(par = uniqueOmegas, fn = stepOmega, indOmegas = omegaRestrictions,
                             objFMM = outMobius, omegaMax = omegaMax,
                             control = list(warn.1d.NelderMead = FALSE))$par
  # warn.1d.NelderMead to suppress Nelder-Mead method warning when used for a single omega.

  # A and M estimates are recalculated by linear regression
  beta <- getBeta(outMobius)
  alpha <- getAlpha(outMobius)
  omega <- uniqueOmegasOptim[omegaRestrictions]

  designMatrix <- matrix(0, ncol = nback, nrow = n)
  for(j in 1:nback){
    designMatrix[,j] <- calculateCosPhi(alpha = alpha[j], beta = beta[j], omega = omega[j], timePoints = timePoints)
  }
  regresion <- lm(vData ~ designMatrix)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]

  # Fitted values
  fittedFMMvalues <- predict(regresion)

  # Residual sum of squares
  SSE <- sum((fittedFMMvalues-vData)^2)

  names(A) <- paste("A", 1:length(A), sep = "")

  nIter <- getNIter(outMobius)

  # Returns an object of class FMM
  outMobiusFinal <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues= fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  )
  return(outMobiusFinal)
}

###############################################################
# Internal function: to fit restricted multicomponent FMM models.
# Nested backfitting algorithm is used for fitting process.
# 'fitFMM_restr' function is used.
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   nback: number of FMM components to be fitted.
#   betaRestrictions: beta's constraint vector.
#   omegaRestrictions: omega's constraint vector.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
#   showProgress: TRUE to display a progress indicator on the console.
# Returns an object of class FMM.
# Note1: alphaGrid and omegaGrid as lists are not supported
# Note2: a previous FMM object refine is not supported
###############################################################
fitFMM_restr_omega_beta<-function(vData, nback, timePoints = seqTimes(length(vData)),
                                  betaRestrictions, omegaRestrictions, maxiter = nback,
                            stopFunction = alwaysFalse, lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                            alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid), omegaMin = 0.0001, omegaMax = 1,
                            omegaGrid = exp(seq(log(omegaMin),log(omegaMax), length.out=lengthOmegaGrid)),
                            numReps = 3, showProgress = TRUE, parallelize = FALSE){
  n <- length(vData)

  if(is.list(alphaGrid)){
    stop("alphaGrid as list not supported when specifying alphaRestrictions")
  }

  if(is.list(omegaGrid)){
    stop("omegaGrid as list not supported when specifying omegaRestrictions")
  }

  # showProgress
  if(showProgress){
    totalMarks <- 50
    partialMarkLength <- 2
    progressHeader<-paste(c("|",rep("-",totalMarks),"|\n|"), collapse = "")
    cat(progressHeader)
    completedPercentage <- 0.00001
    previousPercentage <- completedPercentage
  }

  # omega blocks
  numBlocks <- length(unique(omegaRestrictions))

  # Object initialization
  fittedValuesPerBlock <- matrix(0, ncol = numBlocks, nrow = n)
  fittedFMMPerBlock <- list()
  prevFittedFMMvalues <- NULL

  stopCriteria<-"Stopped by reaching maximum iterations ("
  # Backfitting algorithm: iteration
  for(i in 1:maxiter){

    indBloque <- 1

    # Backfitting algorithm: component for each omega block
    for(j in unique(omegaRestrictions)){
        componentes <- which(omegaRestrictions == j)
        numComponents <- length(componentes)

        # data for block k: difference between vData and all other components fitted values
        # of other blocks
        backFittingData <- vData - apply(as.matrix(fittedValuesPerBlock[,-indBloque]), 1, sum)

        iteraciones <- ifelse(numComponents > 1, min(numComponents + 1, 4), 1)

        # fitting of a block using fitFMM_restr function
        fittedFMMPerBlock[[indBloque]] <- fitFMM_restr(backFittingData, timePoints = timePoints, nback = numComponents,
                                                  betaRestrictions = betaRestrictions[componentes],
                                                  omegaRestrictions = rep(1,numComponents), maxiter = iteraciones,
                                                  lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                                                  alphaGrid = alphaGrid, omegaMax = omegaMax, omegaGrid = omegaGrid,
                                                  numReps = numReps, parallelize = parallelize)
        fittedValuesPerBlock[,indBloque] <- getFittedValues(fittedFMMPerBlock[[indBloque]])

        indBloque <- indBloque + 1

        # showProgress
        if(showProgress){
          completedPercentage <- completedPercentage + 100/(nback*maxiter)
          if(ceiling(previousPercentage) < floor(completedPercentage)){
            progressDone <- paste(rep("=",sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1)
                                             %% partialMarkLength == 0))), collapse = "")
            cat(progressDone)
            previousPercentage <- completedPercentage
          }
        }
    }

    # Check stop criterion
    # Fitted values as sum of all components
    fittedFMMvalues <- apply(fittedValuesPerBlock, 1, sum)

    if(!is.null(prevFittedFMMvalues)){
      if(PV(vData, prevFittedFMMvalues) > PV(vData, fittedFMMvalues)){
        fittedFMMPerBlock <- prevFittedFMMPerBlock
        fittedFMMvalues <- prevFittedFMMvalues
        stopCriteria <- "Stopped by reaching maximum R2 ("
        break
      }
      if(stopFunction(vData, fittedFMMvalues, prevFittedFMMvalues)){
        stopCriteria <- "Stopped by the stopFunction ("
        break
      }
    }
    prevFittedFMMvalues <- fittedFMMvalues
    prevFittedFMMPerBlock <- fittedFMMPerBlock
  }
  nIter <- i

  # showProgress
  if(showProgress){
    if(completedPercentage < 100){
      completedPercentage <- 100
      nMarks <- ifelse(ceiling(previousPercentage) < floor(completedPercentage),
                       sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1)
                            %% partialMarkLength == 0)), 0)
      if (nMarks > 0) {
        cat(paste(rep("=",nMarks), collapse = ""))
        previousPercentage <- completedPercentage
      }
    }
    cat("|\n", stopCriteria, nIter, "iteration(s))", "\n")
  }

  # alpha, beta y omega estimates
  alpha <- unlist(lapply(fittedFMMPerBlock, getAlpha))
  beta <- unlist(lapply(fittedFMMPerBlock, getBeta))
  omega <- unlist(lapply(fittedFMMPerBlock, getOmega))

  # A and M estimates are recalculated by linear regression
  cosPhi <- calculateCosPhi(alpha = alpha, beta = beta, omega = omega, timePoints = timePoints)
  regression <- lm(vData ~ cosPhi)
  M <- coefficients(regression)[1]
  A <- coefficients(regression)[-1]

  # Fitted values
  fittedFMMvalues <- predict(regression)

  # Residual sum of squares
  SSE <- sum((fittedFMMvalues - vData)^2)
  names(A) <- paste("A", 1:length(A), sep = "")

  # Returns an object of class FMM
  outMobius <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  )
  return(outMobius)
}
