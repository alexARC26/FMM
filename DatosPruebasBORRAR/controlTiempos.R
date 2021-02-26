
library(microbenchmark)
library(ggplot2)
library(Rcpp)

data1 <- read.csv("DatosPruebasBORRAR/587098017_Sweep16_Train1_Spike1.csv")
data2 <- read.csv("DatosPruebasBORRAR/595568262_Sweep12_Train1_Spike1.csv")

load("DatosPruebasBORRAR/oldExecTimesFMMunit.Rdata")
load("DatosPruebasBORRAR/newExecTimesFMMunit1.Rdata")
load("DatosPruebasBORRAR/newExecTimesFMMunit2.Rdata")
load("DatosPruebasBORRAR/newExecTimesFMMunit3.Rdata")
load("DatosPruebasBORRAR/newExecTimesFMMunit4.Rdata")

execTimesFMM <- rbind(newExecTimesFMMunit1, newExecTimesFMMunit2, newExecTimesFMMunit3,
                      newExecTimesFMMunit4)


execTimesFMM$method <- factor(execTimesFMM$method, levels = unique(execTimesFMM$method)[c(3,1,4,2)])

plotTimes <- function(execTimesFMM){
  execTimesFMM$time <- execTimesFMM$time/1e09
  ggplot(data = execTimesFMM, aes(y = time, x = expr, fill = method),
         main = "FMM parameters estimation times\n(Non-paralelized, unit) ") +
    geom_violin() +
    facet_wrap(.~method, ncol = 2) + #, scales = "free_x"
    coord_flip()
}
#########################################################################################

library("profvis")

profvis(fitFMM(vData = data1$Voltage))

fitFMM(data1$Voltage)

m2 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 1000, sigmaNoise = 0.05)[[3]]
fitFMM(m2)
m2 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 1000, sigmaNoise = 0.1)[[3]]
fitFMM(m2)
load("m2.Rdata")
m2
m3 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 10000, sigmaNoise = 0.05)[[3]]

#timesNonPar <- microbenchmark(m1 = fitFMM(m1), m2 = fitFMM(m2),
#                              m3 = fitFMM(m3), times = 10)
fitFMM(data1$Voltage)
plot(data1$Voltage)

########################################################################################

fitFMM(data1$Voltage)
newExecTimesFMMunit5 <- microbenchmark(rcppProjectM1 = fitFMM(data1$Voltage, useRcpp = TRUE),
                                       rcppProjectM2 = fitFMM(data2$Voltage), times = 50)
newExecTimesFMMunit6 <- microbenchmark(FMMProjectM1 = fitFMM(data1$Voltage, useRcpp = TRUE),
                                       FMMProjectM2 = fitFMM(data2$Voltage), times = 50)

autoplot(rbind(newExecTimesFMMunit5, newExecTimesFMMunit6)) +
  scale_x_continuous()

data("neuronalAPTrain")
data("neuronalSpike")
generatedFMM1 <- generate_FMM(-20, 10, 2.5, 1.5, 0.1, length.out = 1000)[[3]]
fitFMM(neuronalAPTrain)
fitFMM(neuronalSpike)
plotFMM(fitFMM(neuronalSpike))
