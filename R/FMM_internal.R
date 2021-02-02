###############################################################
# Auxiliary internal functions
# Functions:
#   step1FMM:     M, A and beta initial parameter estimations.
#   bestStep1:    to find the optimal initial parameters estimation.
#   step2FMM:     second step of FMM fitting process.
#   refineFMM:    fitFMM from a previous objectFMM.
#   PV:           percentage of variability explained.
#   PVj:          percentage of variability explained by each component of FMM model.
#   angularmean:  to compute the angular mean.
#   seqTimes:     to build a sequence of equally time points spaced in range [0,2*pi].
###############################################################


###############################################################
# Internal function: to estimate M, A and beta initial parameters
# also returns residual sum of squared (RSS).
# Arguments:
#    alphaOmega: vector of the parameters alpha and omega
#    vData: data to be fitted an FMM model.
#    timePoints: one single period time points.
# Returns a 6-length numerical vector: M, A, alpha, beta, omega and RSS
###############################################################
step1FMM <- function(alphaOmega, vData, timePoints) {

  alpha <- alphaOmega[1]
  omega <- alphaOmega[2]

  parteMobius <- 2*atan(omega*tan((timePoints-alpha)/2))
  t_star <- alpha + parteMobius

  # cosinor model with fixed alpha and omega
  xx<-cos(t_star)
  zz<-sin(t_star)
  fit<-lm((vData)~xx+zz)
  M <-fit$coefficients[1] # intercept
  bb<-fit$coefficients[2] # cos coefficient
  gg<-fit$coefficients[3] # sin coefficient
  phiEst<-atan2(-gg,bb)   # acrophase (phi)
  A<-sqrt(bb^2+gg^2)      # wave amplitude
  beta <- (phiEst+alpha)%%(2*pi)

  dataReg<-cos(beta+parteMobius) # Mobius result without intercept and amplitude

  adj0<-M+A*dataReg # Mobius regression
  RSS<-sum((vData-adj0)^2)/length(timePoints) # residual sum of squares

  devolver <- c(M,A,alpha,beta,omega,RSS)
  return(devolver)

}


###############################################################
# Internal function: to find the optimal initial parameter estimation
# Arguments:
#    vData: data to be fitted an FMM model.
#    step1: a data.frame with estimates of
#           M, A, alpha, beta, omega, RSS as columns.
# Returns the optimal row of step1 argument.
# optimum: minimum RSS with several stability conditions.
###############################################################
bestStep1 <- function(vData,step1){

  # step1 in decreasing order by RSS
  ordenMinRSS <- order(step1[,"RSS"])

  maxVData <- max(vData)
  minVData <- min(vData)
  n <- length(vData)

  # iterative search: go through rows ordered step 1
  #    until the first one that verifies the stability conditions
  condicionContinuar <- TRUE
  i <- 1
  while(condicionContinuar){
    # parameters
    M <- step1[ordenMinRSS[i],"M"]
    A <- step1[ordenMinRSS[i],"A"]
    alpha <- step1[ordenMinRSS[i],"alpha"]
    beta <- step1[ordenMinRSS[i],"beta"]
    omega <- step1[ordenMinRSS[i],"omega"]
    sigma <- sqrt(step1[ordenMinRSS[i],"RSS"]*n/(n-5))

    # stability conditions
    maxi <- M + A
    mini <- M - A
    rest1 <- maxi <= maxVData+1.96*sigma
    rest2 <- mini >= minVData-1.96*sigma

    # it is necessary to check that there are no NA,
    # because it can be an extreme solution
    if(is.na(rest1)) rest1 <- FALSE
    if(is.na(rest2)) rest2 <- FALSE

    if(rest1 & rest2){
      condicionContinuar <- FALSE
    } else {
      i <- i+1
    }

    if(i > nrow(step1)){
      return(NULL)
    }
  }
  return(step1[ordenMinRSS[i],])

}

###############################################################
# Internal function: second step of FMM fitting process
# Arguments:
#   param: M, A, alpha, beta, omega initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omegaMax: max value for omega.
###############################################################
step2FMM <- function(param, vData, timePoints, omegaMax){

  n <- length(timePoints)

  # FMM model
  ffMob<-param[1] + param[2] * cos(param[4]+2*atan2(param[5]*sin((timePoints-param[3])/2),cos((timePoints-param[3])/2)))

  # Residual sum of squares
  RSS<-sum((ffMob - vData)^2)/n
  sigma <- sqrt(RSS*n/(n-5))

  # When amplitude condition is valid, it returns RSS
  # else it returns infinite.
  maxi<-param[1]+param[2]
  mini<-param[1]-param[2]
  rest1 <- maxi <= max(vData)+1.96*sigma
  rest2 <- mini >= min(vData)-1.96*sigma

  # Other integrity conditions that must be met
  rest3 <- param[2] > 0         # A > 0
  rest4 <- param[5] > 0         # omega > 0
  rest5 <- param[5] <= omegaMax # omega <= omegaMax

  if(rest1 & rest2 & rest3 & rest4 & rest5){
    return(RSS)
  }else{
    return(Inf)
    #return(10^10)
  }

}


##################################################################
# Internal function: fitFMM from a FMM object
# Same arguments as fitFMM function and two
# additional arguments:
#   objectFMM: to start fitting from a previous objectFMM
#   staticComponents: vector containing the components that will not be modified
##################################################################
refineFMM <- function(vData, nPeriods = 1, timePoints = NULL,
                      nback = 1, betaRestrictions = 1:nback, omegaRestrictions = 1:nback, maxiter = nback,
                      stopFunction = alwaysFalse, objectFMM = NULL, staticComponents = NULL,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      numReps = 3, showProgress = TRUE, showTime = TRUE, parallelize=FALSE) {

  alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid)
  omegaMax = 1
  omegaGrid = exp(seq(log(0.0001),log(omegaMax),length.out=lengthOmegaGrid))

  if(showTime) time.ini <- Sys.time()

  if(nPeriods > 1){
    n <- length(vData)
    if(n%%nPeriods != 0) stop("Data length is not a multiple of nPeriods")
    M <- matrix(vData,nrow=nPeriods,ncol=n/nPeriods,byrow = TRUE)
    #vDataAnt <- vData
    summarizedData <- apply(M,2,mean)
  } else {
    summarizedData <- vData
  }

  if(is.null(timePoints)){
    timePoints<-seqTimes(length(summarizedData))
  } else {
    if(any(timePoints < 0) | any(timePoints > 2*pi)){
      stop("timePoints must be between 0 and 2*pi")
    }
  }

  if(nback == 1){
    res <- fitFMM_unit(summarizedData, timePoints, lengthAlphaGrid, lengthOmegaGrid, alphaGrid, omegaMax,
                       omegaGrid,numReps)
  } else {
    if(length(unique(betaRestrictions)) == nback & length(unique(omegaRestrictions)) == nback){
      res <- fitFMM_back(summarizedData, timePoints, nback, maxiter, stopFunction, objectFMM, staticComponents,
                         lengthAlphaGrid, lengthOmegaGrid, alphaGrid, omegaMax, omegaGrid, numReps, showProgress)
    } else {
      if(length(unique(omegaRestrictions)) == nback & length(unique(betaRestrictions)) != nback){
        res <- fitFMM_restr_beta(summarizedData, timePoints, nback, betaRestrictions, maxiter, stopFunction, objectFMM,
                                 staticComponents, lengthAlphaGrid, lengthOmegaGrid, alphaGrid, omegaMax, omegaGrid, numReps, showProgress)
      } else {
        if(showProgress)
          warning("showProgress not available when specifying omegaRestrictions.")
        res <- fitFMM_restr(summarizedData, timePoints, nback, betaRestrictions, omegaRestrictions, maxiter, stopFunction, objectFMM,
                            staticComponents, lengthAlphaGrid, lengthOmegaGrid, alphaGrid, omegaMax, omegaGrid, numReps, parallelize)
      }
    }
  }

  if(showTime){
    time.end <- Sys.time()
    cat(time.end-time.ini)
  }

  res@nPeriods <- nPeriods
  res@data <- vData

  ordenVariabilidad <- order(res@R2,decreasing = TRUE)

  res@A <- res@A[ordenVariabilidad]
  res@alpha <- res@alpha[ordenVariabilidad]
  res@beta <- res@beta[ordenVariabilidad]
  res@omega <- res@omega[ordenVariabilidad]
  res@R2 <- res@R2[ordenVariabilidad]

  return(res)

}

##################################################################
# Internal function: to calculate the percentage of variability explained by the FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values.
##################################################################
PV <- function(vData,pred){
  meanVData <- mean(vData)
  return(1 - sum((vData-pred)^2)/sum((vData-meanVData)^2))
}

##################################################################
# Internal function: to calculate the percentage of variability explained by each component of FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   alpha, beta, omega: vectors of corresponding parameter estimates.
##################################################################
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



##################################################################
# Internal function: to build a sequence of equally time points spaced
#                    in range [0,2*pi].
# Arguments:
#   n: secuence length.
##################################################################
seqTimes <- function(n){
  timePoints<-seq(0,2*pi,by=2*pi/n)
  timePoints<-timePoints[-length(timePoints)]
  return(timePoints)
}


##################################################################
# Internal function: to compute the angular mean.
# Arguments:
#   angles: input vector of angles.
##################################################################
#Computes the Angular Mean
angularmean <- function(angles){
  n <- length(angles)
  a.mean <- atan2(sum(sin(angles)),sum(cos(angles)))
  return(a.mean)
}
