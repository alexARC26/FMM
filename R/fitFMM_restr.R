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
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
#   parallelize: TRUE to use parallelized procedure to fit restricted FMM model.
# Returns an object of class FMM.
###############################################################
fitFMM_restr<-function(vData, timePoints = seqTimes(length(vData)), nback,
                      betaRestrictions, omegaRestrictions, maxiter = nback,
                      stopFunction = alwaysFalse, objectFMM = NULL, staticComponents = NULL,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMax = 1,
                      omegaGrid = exp(seq(log(0.001),log(omegaMax), length.out = lengthOmegaGrid)),
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

  if(parallelize){
    nCores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(nCores, outfile="")
    doParallel::registerDoParallel(cl)
  }

  # External loop on the omega grid, setting its value
  objectFMMList <- foreach::foreach(omegas = iterators::iter(omegasIter, by="row")) %dopar% {

    omegas <- as.numeric(omegas)

    # Object initialization
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
          if(j != k){
            vDataAjuste <- vDataAjuste - predichosComponente[[k]]
          }
        }

        # component j fitting using fitFMM_unit_restr function
        ajusteComponente[[j]] <- fitFMM_unit_restr(vDataAjuste,omegas[j], timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                             alphaGrid = alphaGrid[[j]], numReps = numReps)
        predichosComponente[[j]] <- getFittedValues(ajusteComponente[[j]])

      }

      # Check stop criterion
      # Fitted values as sum of all components
      adjMob <- rep(0,n)
      for(j in 1:nback){
        adjMob <- adjMob + predichosComponente[[j]]
      }
      if(!is.null(prevAdjMob)){

        if(PV(vData,prevAdjMob) > PV(vData,adjMob)){
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

    nIter <- i

    # alpha, beta y omega estimates
    alpha <- rep(0,nback)
    beta <- rep(0,nback)
    omega <- rep(0,nback)
    for(j in 1:nback){
        alpha[j] <- getAlpha(ajusteComponente[[j]])
        beta[j] <- getBeta(ajusteComponente[[j]])
        omega[j] <- getOmega(ajusteComponente[[j]])
    }

    ## Check if the solution is valid
    # beta restrictions: calculate angular mean of beta parameters
    # the nearest betas are chosen
    RestBeta<-beta
    elegidos <- rep(0,length(betaRestrictions))
    vCompleto <- 1:length(betaRestrictions)
    for(indRes in unique(betaRestrictions)){
      numComponents <- sum(betaRestrictions == indRes)
      primero <- vCompleto[elegidos == 0][1]
      distanciaAbs <- abs(RestBeta-RestBeta[primero])
      distanciaAbs[elegidos == 1] <- Inf
      implicados <- order(distanciaAbs)[1:numComponents]
      elegidos[implicados] <- 1
      RestBeta[implicados] <- angularmean(RestBeta[implicados])%%(2*pi)
    }

    # A and M estimates are recalculated by linear regression
    cos.phi <- list()
    for(j in 1:nback){
      cos.phi[[j]] <- cos(RestBeta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
    }
    M <- matrix(unlist(cos.phi),ncol=nback)
    regresion <- lm(vData ~ M)
    M <- coefficients(regresion)[1]
    A <- coefficients(regresion)[-1]

    SSE.k <- 0
    if(sum(A < 0,na.rm=TRUE) > 0){
      SSE.k <- 10^6
    }


    # Fitted values:
    cos.phi <- list()
    for(j in 1:nback){
      cos.phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
    }
    M <- matrix(unlist(cos.phi),ncol=nback)
    regresion <- lm(vData ~ M)
    M <- coefficients(regresion)[1]
    A <- coefficients(regresion)[-1]

    adjMob <- predict(regresion)

    # Residual sum of squares
    SSE <- sum((adjMob-vData)^2) + SSE.k

    names(A) <- paste("A", 1:length(A), sep="")

    if(parallelize) parallel::stopCluster(cl)

    # Returns an object of class FMM
    outMobius <- FMM(
      M = M,
      A = A,
      alpha = alpha,
      beta = beta,
      omega = omega,
      timePoints = timePoints,
      summarizedData = vData,
      fittedValues = adjMob,
      SSE = SSE,
      R2 = PVj(vData, timePoints, alpha, beta, omega),
      nIter = nIter
    )

    return(outMobius)
  }

  # We keep the best solution
  SSEs <- rep(NA,length(objectFMMList))
  for(i in 1:length(SSEs)){
    SSEs[i] <- getSSE(objectFMMList[[i]])
  }
  imin <- which.min(SSEs)
  outMobius <- objectFMMList[[imin]]

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

  elegidos <- rep(0,length(betaRestrictions))
  vCompleto <- 1:length(betaRestrictions)
  for(indRes in unique(betaRestrictions)){
    numComponents <- sum(betaRestrictions == indRes)
    primero <- vCompleto[elegidos == 0][1]
    distanciaAbs <- abs(beta-beta[primero])
    distanciaAbs[elegidos == 1] <- Inf
    implicados <- order(distanciaAbs)[1:numComponents]
    elegidos[implicados] <- 1
    beta[implicados] <- angularmean(beta[implicados])%%(2*pi)
  }

  cos.phi <- list()
  for(j in 1:nback){
    cos.phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  M <- matrix(unlist(cos.phi),ncol=nback)
  regresion <- lm(vData ~ M)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]

  # Fitted values
  adjMob <- predict(regresion)

  # Residual sum of squares
  SSE <- sum((adjMob-vData)^2)

  names(A) <- paste("A",1:length(A),sep="")

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
    fittedValues= adjMob,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  )

  return(outMobiusFinal)

}

###############################################################
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
###############################################################
stepOmega <- function(uniqueOmegas, indOmegas, objFMM, omegaMax){

  timePoints <- getTimePoints(objFMM)
  n <- length(timePoints)
  M <- getM(objFMM)
  A <- getA(objFMM)
  alpha <- getAlpha(objFMM)
  beta <- getBeta(objFMM)
  omega <- uniqueOmegas[indOmegas]
  m <- length(alpha)
  vData <- getSummarizedData(objFMM)

  # FMM fitting
  wave <- list()
  for(j in 1:m){
    wave[[j]] <- A[j]*cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  ffMob<-rep(M,n)
  for(j in 1:length(wave)){
    ffMob <- ffMob+wave[[j]]
  }

  # Residual sum of squares
  RSS<-sum((ffMob - vData)^2)/n
  sigma <- sqrt(RSS*n/(n-5))

  # Other integrity conditions that must be met
  rest1 <- all(uniqueOmegas <= omegaMax)
  rest2 <- all(uniqueOmegas >= 0)
  if(rest1 & rest2){
    return(RSS)
  }else{
    return(Inf)
  }

}


###############################################################
# Internal function: to fit monocomponent FMM models with fixed omega.
# Arguments:
#   vData: data to be fitted an FMM model.
#   omega: value of the omega parameter.
#   timePoints: one single period time points.
#   lengthAlphaGrid: precision of the grid of alpha parameter.
#   alphaGrid: grid of alpha parameter.
#   numReps: number of times the alpha-omega grid search is repeated.
# Returns an object of class FMM.
###############################################################
fitFMM_unit_restr<-function(vData, omega, timePoints = seqTimes(length(vData)),
                            lengthAlphaGrid = 48, alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                            numReps = 3){

  n <- length(vData)

  ## Step 1: initial values of M, A, alpha, beta and omega
  # alpha and omega are fixed and cosinor model is used to calculate the rest of the parameters.
  # step1FMM function is used to make this estimate
  grid <- expand.grid(alphaGrid, omega)
  step1 <- t(apply(grid, 1, step1FMM, vData = vData, timePoints = timePoints))
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
    nelderMead <- optim(par = bestPar[1:4], fn = step2FMM_restr, vData = vData, timePoints = timePoints, omega = omega)
    parFinal <- c(nelderMead$par, omega)
  } else {
    parFinal <- c(bestPar[1:4],omega)
  }

  names(parFinal) <- c("M","A","alpha","beta","omega")
  adjMob <- parFinal["M"] + parFinal["A"]*cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((adjMob-vData)^2)

  # alpha and beta between 0 and 2pi
  parFinal[3] <- parFinal[3]%%(2*pi)
  parFinal[4] <- parFinal[4]%%(2*pi)

  # the grid search is repeated numReps
  numReps <- numReps - 1
  while(numReps > 0){

    # new grid for alpha between 0 and 2pi
    nAlphaGrid <- length(alphaGrid)
    amplitudeAlphaGrid <- 1.5*mean(diff(alphaGrid))
    alphaGrid <- seq(parFinal[3]-amplitudeAlphaGrid,parFinal[3]+amplitudeAlphaGrid,length.out = nAlphaGrid)
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
      nelderMead <- optim(par = bestPar[1:4], fn = step2FMM_restr, vData = vData, timePoints = timePoints, omega = omega)
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
  adjMob <- parFinal["M"] + parFinal["A"]*cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
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
    R2 = 0,
    nIter = 0
  )

  return(outMobius)
}


###############################################################
# Internal function: second step of FMM fitting process with fixed omega
# Arguments:
#   param: M, A, alpha, beta initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omega: fixed value of omega.
###############################################################
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
fitFMM_restr_omega_beta<-function(vData, timePoints = seqTimes(length(vData)), nback,
                            betaRestrictions, omegaRestrictions, maxiter = nback,
                            stopFunction = alwaysFalse, lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                            alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid), omegaMin = 0.001, omegaMax = 1,
                            omegaGrid = exp(seq(log(omegaMin),log(omegaMax),length.out=lengthOmegaGrid)),
                            numReps = 3, showProgress = TRUE){

  n <- length(vData)

  if(is.list(alphaGrid)){
    stop("alphaGrid as list not supported when specifying alphaRestrictions")
  }

  if(is.list(omegaGrid)){
    stop("omegaGrid as list not supported when specifying omegaRestrictions")
  }

  # showProgress
  if(showProgress){
    marcasTotales <- 50
    granoInforme <- 2
    cat("|")
    for(m in 1:marcasTotales) cat("-")
    cat("|\n")
    cat("|")
    porcentajeCompletado <- 0.00001
    porcentajeAntes <- porcentajeCompletado
  }

  # omega blocks
  numBlocks <- length(unique(omegaRestrictions))

  # Object initialization
  predichosBloque <- list()
  ajusteBloque <- list()
  prevAdjMob <- NULL

  for(i in 1:numBlocks){
    predichosBloque[[i]] <- rep(0,n)
  }

  # Backfitting algorithm: iteration
  for(i in 1:maxiter){

    indBloque <- 1

    # Backfitting algorithm: component for each omega block
    for(j in unique(omegaRestrictions)){

        componentes <- which(omegaRestrictions == j)
        numComponents <- length(componentes)

        # data for block k: difference between vData and all other components fitted values
        # of other blocks
        vDataAjuste <- vData
        for(k in 1:numBlocks){
          if(k != indBloque){
            vDataAjuste <- vDataAjuste - predichosBloque[[k]]
          }
        }

        if(numComponents > 1){
          iteraciones <- min(numComponents + 1, 4)
        } else {
          iteraciones <- 1
        }

        # fitting of a block using fitFMM_restr function
        ajusteBloque[[indBloque]] <- fitFMM_restr(vDataAjuste, timePoints = timePoints, nback = numComponents,
                                                  betaRestrictions = betaRestrictions[componentes],
                                                  omegaRestrictions = rep(1,numComponents), maxiter = iteraciones,
                                                  lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                                                  alphaGrid = alphaGrid, omegaMax = omegaMax, omegaGrid = omegaGrid,
                                                  numReps = numReps)

        predichosBloque[[indBloque]] <- getFittedValues(ajusteBloque[[indBloque]])

        indBloque <- indBloque + 1

        # showProgress
        if(showProgress){
          porcentajeCompletado <- porcentajeCompletado + (100*numComponents)/(nback*maxiter)
          if(ceiling(porcentajeAntes) < floor(porcentajeCompletado)){
            numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
          } else {
            numMarcas <- 0
          }
          if (numMarcas > 0) {
            for(m in 1:numMarcas) cat("=")
            porcentajeAntes <- porcentajeCompletado
          }
        }
    }

    # Check stop criterion
    # Fitted values as sum of all components
    adjMob <- rep(0,n)
    for(j in 1:numBlocks){
      adjMob <- adjMob + predichosBloque[[j]]
    }
    if(!is.null(prevAdjMob)){

      if(PV(vData,prevAdjMob) > PV(vData,adjMob)){
        ajusteBloque <- ajusteBloqueAnt
        adjMob <- prevAdjMob
        break
      }

      if(stopFunction(vData,adjMob,prevAdjMob)){
        break
      }

    }

    prevAdjMob <- adjMob
    ajusteBloqueAnt <- ajusteBloque

  }
  nIter <- i

  # showProgress
  if(showProgress){
    if(porcentajeCompletado < 100){
      porcentajeCompletado <- 100
      if(ceiling(porcentajeAntes) < floor(porcentajeCompletado)){
        numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
      } else {
        numMarcas <- 0
      }
      if (numMarcas > 0) {
        for(m in 1:numMarcas) cat("=")
        porcentajeAntes <- porcentajeCompletado
      }
    }
    cat("|\n")
    if(i == maxiter){
      if(i == 1){
        cat("Stopped by reaching maximum iterations (",i ,"iteration )","\n")
      } else {
        cat("Stopped by reaching maximum iterations (",i ,"iterations )","\n")
      }
    } else {
      if(i == 1){
        cat("Stopped by the stopFunction (",i ,"iteration )","\n")
      } else {
        cat("Stopped by the stopFunction (",i ,"iterations )","\n")
      }
    }
  }

  # alpha, beta y omega estimates
  alpha <- c()
  beta <- c()
  omega <- c()
  for(j in 1:numBlocks){
    alpha <- c(alpha,getAlpha(ajusteBloque[[j]]))
    beta <- c(beta,getBeta(ajusteBloque[[j]]))
    omega <- c(omega,getOmega(ajusteBloque[[j]]))
  }

  # beta restrictions: calculate angular mean of beta parameters
  # the nearest betas are chosen
  elegidos <- rep(0,length(betaRestrictions))
  vCompleto <- 1:length(betaRestrictions)
  for(indRes in unique(betaRestrictions)){
    numComponents <- sum(betaRestrictions == indRes)
    primero <- vCompleto[elegidos == 0][1]
    distanciaAbs <- abs(beta-beta[primero])
    distanciaAbs[elegidos == 1] <- Inf
    implicados <- order(distanciaAbs)[1:numComponents]
    elegidos[implicados] <- 1
    beta[implicados] <- angularmean(beta[implicados])%%(2*pi)
  }

  # A and M estimates are recalculated by linear regression
  cos.phi <- list()
  for(j in 1:nback){
    cos.phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  M <- matrix(unlist(cos.phi),ncol=nback)
  regresion <- lm(vData ~ M)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]

  # Fitted values
  adjMob <- predict(regresion)

  # Residual sum of squares
  SSE <- sum((adjMob-vData)^2)

  names(A) <- paste("A", 1:length(A),sep="")

  # Returns an object of class FMM
  outMobius <- FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = adjMob,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  )

  return(outMobius)
}
