# Extract individual contribution to the fitted values of each FMM wave
#
# Arguments:
#   objFMM: object of class FMM.
# Returns a list with as many elements as there are components in the model.
extractWaves <- function(objFMM){
   nComponents <- length(getAlpha(objFMM))
   timePoints <- getTimePoints(objFMM)
   firstValue<-getData(objFMM)[1]

   predicted <- list()

   for(i in 1:nComponents){
     predictedi <- getA(objFMM)[i]*cos(getBeta(objFMM)[i] + 2*atan(getOmega(objFMM)[i]*tan((timePoints-getAlpha(objFMM)[i])/2)))
     predicted[[i]] <- predictedi - predictedi[1] + firstValue
   }

   return(predicted)
}
