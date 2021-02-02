# Plot fitted FMM models
#
# Arguments:
#   objFMM: object of class FMM.
#   components: TRUE to plot each component separately.
#   plotAlongPeriods: TRUE to plot more than 1 period.
#   use_ggplot2: TRUE to use ggplot2 package.
#   legendInComponentsPlot: TRUE to indicate if the a legend should be plotted in the component plot.
#   textExtra: Extra text to be added to the titles of the plot.
plotFMM <- function(objFMM, components=FALSE, plotAlongPeriods=FALSE,
                    use_ggplot2=FALSE, legendInComponentsPlot = TRUE, textExtra = ""){

  nPeriods=getNPeriods(objFMM)
  if(nPeriods>1){
    if(plotAlongPeriods & !components){
      vData <- getData(objFMM)
    }else{
      vData <- getSummarizedData(objFMM)
    }
  }else{vData <- getData(objFMM)}
  nObs<-length(vData)

  if(plotAlongPeriods & !components){
    timePoints <- getTimePoints(objFMM)
    timePoints <- rep(timePoints,nPeriods)
  }else{
    timePoints <- getTimePoints(objFMM)
  }

  col = 2:(length(getAlpha(objFMM))+1)


  significantTimePoints<-round(c(1,nObs*0.25, nObs*0.5, nObs*0.75, nObs))

  # Components plot: if there is more than one period, just the data from the first period will be plotted
  if(components){
    title <- ifelse(textExtra != "", paste("Components FMM",textExtra,sep = " - "),"Components FMM")
    nComponents <- length(getAlpha(objFMM))
    componentNames<-paste("Wave ",1:nComponents, sep="")

    predicted<-extractWaves(objFMM)

    minValue <- min(sapply(predicted,min))
    maxValue <- max(sapply(predicted,max))

    if(!use_ggplot2){
      plot(1:nObs,vData,ylim=c(minValue,maxValue),xlab="Time",ylab="Response",main=title,type="n",xaxt = "n")
      for(i in 1:nComponents){
        points(1:nObs,predicted[[i]],type="l",lwd=2,col=col[i])
      }
      axis(1, las = 1,
           at = significantTimePoints,
           labels = parse(text=paste("t[",significantTimePoints, "]", sep = "")))
      if(legendInComponentsPlot) legend("topright",legend=componentNames,col=col,lty=1)
    } else {
      requireNamespace("ggplot2", quietly = TRUE)
      requireNamespace("RColorBrewer", quietly = TRUE)

      df<-data.frame("Time"=rep(1:length(timePoints),nComponents),
                     "Response"=unlist(predicted),
                     "Components"=rep(componentNames,each=nObs))

      # With more than 9 components, the selection of colors must be expanded
      if(nComponents>9){
        colorsForComponents <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nComponents)
      }else{
        colorsForComponents <- ifelse(rep(nComponents>3,nComponents),RColorBrewer::brewer.pal(nComponents, "Set1"),RColorBrewer::brewer.pal(3, "Set1"))
      }

      plot<-ggplot2::ggplot(data=df, ggplot2::aes_(x=~Time, y=~Response,
                                                  group=~Components, color=~Components)) +
        ggplot2::geom_line(ggplot2::aes_(color=~Components),
                  size=1.3,lineend = "round",linejoin = "round")+
        ggplot2::scale_color_manual(values=colorsForComponents)+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = ifelse(legendInComponentsPlot,"bottom","none"))+
        ggplot2::labs(title = title)+
        ggplot2::scale_x_continuous(
          breaks = significantTimePoints,
          labels = function(x) parse(text=paste("t[",x,"]")))
      return(plot)

    }

  } else {
    title <- ifelse(textExtra != "", paste("Fitted FMM model",textExtra,sep = " - "),"Fitted FMM model")

    if(!use_ggplot2){
      plot(1:nObs,vData,xlab="Time",ylab="Response",main=title,xaxt = "n")
      if(plotAlongPeriods){
        points(1:nObs,rep(getFittedValues(objFMM),nPeriods),type="l",col=2,lwd=2)
      }else{
        points(1:nObs,getFittedValues(objFMM),type="l",col=2,lwd=2)
      }
      axis(1, las = 1,
           at = significantTimePoints,
           labels = parse(text=paste("t[",significantTimePoints, "]", sep = "")))
    } else {
      requireNamespace("ggplot2", quietly = TRUE)

      if(plotAlongPeriods){
        adjustedModel<-rep(getFittedValues(objFMM),nPeriods)
      }else{
        adjustedModel<-getFittedValues(objFMM)
      }

      fittedData<-data.frame("Time"=1:nObs,
                             "fitted_FMM"=adjustedModel,
                             "Response"=vData)

      plot<-ggplot2::ggplot(data=fittedData,
                            ggplot2::aes_(x=~Time, y=~Response, color=1)) +
        ggplot2::geom_point(size=2,color="grey65", shape=21, stroke=1.1)+
        ggplot2::geom_path(ggplot2::aes_(x=~Time, y=~fitted_FMM, color="FMM", position=NULL),
                  size=2,lineend = "round",linejoin = "round")+
        ggplot2::labs(title = title)+
        ggplot2::scale_color_manual(values="red")+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "none")+
        ggplot2::scale_x_continuous(
          breaks = significantTimePoints,
          labels = function(x) parse(text=paste("t[",x,"]")))

      return(plot)
    }

  }

}
