# Fit FMM model
#
# Arguments:
#     vData: A numeric vector which contains the data to be fitted an FMM model.
#     nPeriods: A numeric value specifying the number of periods at which
#               \code{vData} is observed.
#     timePoints: A numeric vector containing the time points at which each data
#                 of one single period is observed.
#                 The default value is NULL, in which case they are equally
#                 spaced in range [0,2*pi].
#                 It must be between 0 to 2*pi.
#     nback: Number of FMM components to be fitted. Its default value is 1.
#     betaRestrictions: An integer vector of length nback indicating which FMM
#                       waves are constrained
#                       to have equal beta parameters.
#     omegaRestrictions: An integer vector of length \code{nback} indicating
#                        which FMM waves are constrained
#                        to have equal omega parameters.
#     maxiter: Maximum number of iterations for the backfitting algorithm.
#              By default, it is setting at nback.
#     stopFunctions: Function to check the criterion convergence for the
#                    backfitting algorithm.
#     lengthAlphaGrid: Precission of the grid of alpha in the search of the best
#                      model. By default it's established at 48 possible values
#                      of alpha, equally spaced between 0 and 2*pi.
#     lengthOmegaGrid: Precission of the grid of omega in the search of the best
#                      model. By default it's established at 24 possible values
#                      of omega, equally spaced between 0 and 1 in a logarithmic
#                      way.
#     numReps: Number of times the fitting is repeated.
#     showProgress: TRUE to display a progress indicator on the console.
#     showTime: TRUE to display execution time on the console.
#     parallelize: TRUE to use parallelized procedure to fit restricted
#                  FMM model.
#     useRcpp: TRUE to use Rcpp function.
fitFMM <- function(vData, nPeriods = 1, timePoints = NULL,
                   nback = 1, betaRestrictions = 1:nback,
                   omegaRestrictions = 1:nback, maxiter = nback,
                   stopFunction = alwaysFalse,
                   lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                   numReps = 3, showProgress = TRUE, showTime = TRUE,
                   parallelize = FALSE, useRcpp = FALSE){

  alphaGrid <- seq(0,2*pi,length.out = lengthAlphaGrid)
  omegaMax <- 1
  omegaGrid <- exp(seq(log(0.0001),log(omegaMax),length.out = lengthOmegaGrid))
  staticComponents <- NULL
  objectFMM <- NULL

  betaRestrictions <- sort(betaRestrictions)
  omegaRestrictions <- sort(omegaRestrictions)

  if(showTime) time.ini <- Sys.time()

  # If data has more than one period, it must be summarized
  if(nPeriods > 1){
    n <- length(vData)
    if(n %% nPeriods != 0) stop("Data length is not a multiple of nPeriods")
    dataMatrix <- matrix(vData, nrow = nPeriods, ncol = n/nPeriods, byrow = TRUE)
    #vDataAnt <- vData
    summarizedData <- apply(dataMatrix,2,mean)
  } else {
    summarizedData <- vData
  }

  # Generation of the time points
  if(is.null(timePoints)){
    timePoints<-seqTimes(length(summarizedData))
  } else {
    if(any(timePoints < 0) | any(timePoints > 2*pi)){
      stop("timePoints must be between 0 and 2*pi")
    }
    if(length(timePoints) != length(summarizedData)){
      stop("timePoints must have the same length as one-period data")
    }
  }

  # used apply function for compute FMM models
  usedApply_Cluster <- getApply(parallelize)
  usedApply <- usedApply_Cluster[[1]]

  if(nback == 1){
    fittedFMM <- fitFMM_unit(vData = summarizedData, timePoints = timePoints,
                       lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                       alphaGrid = alphaGrid, omegaMax = omegaMax, omegaGrid = omegaGrid,
                       numReps = numReps, usedApply = usedApply, useRcpp)

  } else {
    if(length(unique(betaRestrictions)) == nback &
       length(unique(omegaRestrictions)) == nback){
      fittedFMM <- fitFMM_back(summarizedData,timePoints, nback, maxiter,stopFunction,
                         objectFMM, staticComponents, lengthAlphaGrid,
                         lengthOmegaGrid, alphaGrid, omegaMax, omegaGrid,
                         numReps, showProgress, usedApply = usedApply, useRcpp )
    } else {
      if(length(unique(omegaRestrictions)) == nback &
         length(unique(betaRestrictions)) != nback){
        fittedFMM <- fitFMM_restr_beta(summarizedData, timePoints, nback,
                                 betaRestrictions, maxiter, stopFunction,
                                 objectFMM, staticComponents, lengthAlphaGrid,
                                 lengthOmegaGrid, alphaGrid, omegaMax,omegaGrid,
                                 numReps, showProgress)
      } else {
        fittedFMM <- fitFMM_restr_omega_beta(vData,timePoints, nback,betaRestrictions,
                                       omegaRestrictions, maxiter, stopFunction,
                                       lengthAlphaGrid, lengthOmegaGrid,
                                       alphaGrid, omegaMax, omegaGrid,numReps,
                                       showProgress)
      }
    }
  }

  cluster <- usedApply_Cluster[[2]]
  if(!is.null(cluster)) parallel::stopCluster(cluster)


  if(showTime & showProgress){
    time.end <- Sys.time()
    cat(time.end-time.ini)
  }

  fittedFMM@nPeriods <- nPeriods
  fittedFMM@data <- vData

  explainedVarOrder <- order(fittedFMM@R2,decreasing = TRUE)

  fittedFMM@A <- fittedFMM@A[explainedVarOrder]
  fittedFMM@alpha <- fittedFMM@alpha[explainedVarOrder]
  fittedFMM@beta <- fittedFMM@beta[explainedVarOrder]
  fittedFMM@omega <- fittedFMM@omega[explainedVarOrder]
  fittedFMM@R2 <- fittedFMM@R2[explainedVarOrder]

  # Restricted algorithm may find models with A<0
  needFix <- which(fittedFMM@A < 0)
  if(length(needFix)>0) {
    stop("Invalid solution: check function input parameters.")
  }
  return(fittedFMM)
}




