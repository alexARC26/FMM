\docType{data}
\name{neuronalAPTrain}
\alias{neuronalAPTrain}
\title{Neuronal AP Train Data simulated with Hodgkin-Huxley model}
\format{
	A numeric vector.
}
\source{
	NeuroDynex Documentation, <https://lcn-neurodynex-exercises.readthedocs.io/en/latest/#>
}
\usage{
data(neuronalSpike)
}
\description{
	Voltage data in mV simulated with Hodgkin Huxley model (parameters: C=1,
	gNa=64, gK=39, gL=0.12, vK=-12, vNa=115, vL=10.6, bar(alphaN)=0.85,
	bar(betaN)=1.3, bar(alphaM)=1, bar(betaM)=1.3, bar(alphaH)=1,
	bar(betaH)=1) and applied current of 4.5 microA 1 millisecond. The
	simulation has been done with a modified 'NeuroDynex Python' module.
}
\examples{
data(neuronalAPTrain)
str(neuronalAPTrain)
}
\references{
	Wulfram Gerstner, Werner M. Kistler, Richard Naud,
	and Liam Paninski (2014). Neuronal Dynamics: From Single Neurons to
	Networks and Models of Cognition.
	([Online Book](https://neuronaldynamics.epfl.ch/online/))
}
\keyword{datasets}
