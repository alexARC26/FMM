

################################################################################
# Log-likelihood function for a single wave
# Arguments:
#   parameters:  alpha, omega 
#   vData: data to be fitted an FMM model.
#   varErr: estimated variance of the error (variance of residuals when omega 
#           and alpha are ML).
################################################################################

logLikFMMWave <- function(parameters, data, varErr){
  alpha <- parameters[1]%%(2*pi)
  omega <- parameters[2]
  n <- length(data)
  t <- seq(0, 2*pi, length.out = n)
  t2 <- alpha + 2*atan(omega*tan((t-alpha)/2))
  dM <- cbind(rep(1, n), cos(t2), sin(t2))
  mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%data
  M <- mean(data) - mDeltaGamma[2]*mean(cos(t2)) - mDeltaGamma[3]*mean(sin(t2))
  yFit <- M + mDeltaGamma[2]*cos(t2) + mDeltaGamma[3]*sin(t2)
  return(unname(-sum((data-yFit)^2)/(2*varErr)))
}




























