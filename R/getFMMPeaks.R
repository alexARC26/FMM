# Estimate peak and trough time to a fitted FMM model
#
# Arguments:
#   objFMM: object of class FMM.
#   timePointsIn2pi: TRUE to return peak and trough times is the interval [0,2*pi]
# Returns a list with as many elements as there are components in the model.

#####################################################
##### Function to  #####
## Input argument: FMM object
getFMMPeaks <- function(objFMM, timePointsIn2pi = TRUE) {

  M <- getM(objFMM)
  A <- getA(objFMM)
  alpha <- getAlpha(objFMM)
  beta <- getBeta(objFMM)
  omega <- getOmega(objFMM)

  data<- getData(objFMM)
  if(getNPeriods(objFMM) > 1) data<-getSummarizedData(objFMM)

  nData <- length(data)


  nComp <- length(alpha)

  # timePoints estimation
  peakU<-(alpha+2*atan2(1/omega*sin(-beta/2),cos(-beta/2)))%%(2*pi)
  peakL<-(alpha+2*atan2(1/omega*sin((pi-beta)/2),cos((pi-beta)/2)))%%(2*pi)

  if(timePointsIn2pi){
    tpeakU<-peakU
    tpeakL<-peakL
  }else{
    tpeakU<-peakU*nData/(2*pi) + 1
    tpeakL<-peakL*nData/(2*pi) + 1
  }

  # signal estimation
  phU <- lapply(1:nComp,function(k) (peakU[k]-alpha)/2)
  phL <- lapply(1:nComp,function(k) (peakL[k]-alpha)/2)
  ZU <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phU[[k]])))))
  ZL <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phL[[k]])))))
  names(ZU)<-names(ZL)<-NULL

  return(list(tpeakU=tpeakU, tpeakL=tpeakL, ZU=ZU, ZL=ZL))
}
