###############################################################
# Internal functions to check the criterion convergence for the backfitting algorithm.
# Functions:
#   alwaysFalse:  returns FALSE to force maxiter iterations.
#   R2: check the converge criterion based on explained variability.
###############################################################

###############################################################
# Internal function: FALSE
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values from current iteration.
#   prevPred: fitted values from previous iteration
# Returns FALSE to force maxiter iterations.
###############################################################
alwaysFalse <- function(vData,pred,prevPred){
  return(FALSE)
}

###############################################################
# Internal function: to check if the convergence criterion
#                    based on the difference between the
#                    explained variability in two consecutive
#                    iterations, is reached.
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values from current iteration.
#   prevPred: fitted values from previous iteration.
#   difMax: smallest observable difference between iterations
#           to stop backfitting algorithm.
# Returns TRUE when the converge criterion is reached.
###############################################################
R2 <- function(vData,pred,prevPred,difMax = 0.001){
  funDevolver <- function(vData,pred,prevPred){
    R2.prev <- PV(vData,prevPred)
    R2.now <- PV(vData,pred)
    R2.dif <- R2.now - R2.prev

    if(R2.dif < difMax) return(TRUE)
    return(FALSE)
  }
  return(funDevolver)
}

