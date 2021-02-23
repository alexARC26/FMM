
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
execTimesFMM$time <- execTimesFMM$time/1e09

execTimesFMM$method <- factor(execTimesFMM$method, levels = unique(execTimesFMM$method)[c(3,1,4,2)])
ggplot(data = execTimesFMM, aes(y = time, x = expr, fill = method),
       main = "FMM parameters estimation times\n(Non-paralelized, unit) ") +
  geom_violin() +
  facet_wrap(.~method, ncol = 2) + #, scales = "free_x"
  coord_flip()

#########################################################################################

# OPTIMIZACIÓN DE FMM_unit

library("profvis")

profvis(fitFMM(vData = data1$Voltage))

fitFMM(data2$Voltage)

m1 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 500, sigmaNoise = 0.05)[[3]]
m2 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 1000, sigmaNoise = 0.05)[[3]]
m3 <- generate_FMM(-3, 15, 2.5, 0.3, 0.04, length.out = 10000, sigmaNoise = 0.05)[[3]]

#timesNonPar <- microbenchmark(m1 = fitFMM(m1), m2 = fitFMM(m2),
#                              m3 = fitFMM(m3), times = 10)

fitFMM(m2)
fitFMM(data1$Voltage)
plot(data1$Voltage)


