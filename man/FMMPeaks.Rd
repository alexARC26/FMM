\name{FMMPeaks}
\alias{FMMPeaks}
\title{Peak and trough times and signal values}
\usage{
FMMPeaks(objFMM, timePointsIn2pi = TRUE)
}
\description{
  \code{FMMPeaks()} can be used to estimate interesting parameters of the FMM model specially useful in applications. The function estimates peak and trough times and signal values at those times for each component of the model.
}
\arguments{
	\item{objFMM}{Object of class \code{'FMM'}.}
	\item{timePointsIn2pi}{\code{TRUE} to return peak and trough times in the [0,2*pi] interval.
	                       When \code{timePointsIn2pi = FALSE} the positions of peak and trough times are returned.
	                       Its default value is \code{TRUE}.}
}

\value{
A list with the following components is returned:
  \item{tpeakU}{a numeric vector with the time points at which the peak of each wave is estimated.}
  \item{tpeakL}{a numeric vector with the time points at which the trough of each wave is estimated.}
  \item{ZU}{a numeric vector with the estimated signal peak values of each wave.}
  \item{ZL}{a numeric vector with the estimated signal trough values of each wave.}
}

\references{
Rueda C, Larriba Y, Peddada SD (2019).
   Frequency Modulated Moebius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences.
   \emph{Scientific reports}, \bold{9} (1), 18701. \url{https://www.nature.com/articles/s41598-019-54569-1}
}

\examples{
## Generate example data:
fmm2.data <- generate_FMM(0, 2, c(1.5, 3.4), c(0.2, 2.3), c(0.1, 0.2),
                    plot = FALSE, outvalues = TRUE,
                    sigmaNoise = 0.5) # add a gaussian noise with sigma = 0.5

## Fit the FMM model with nback = 2 components
## fit is an object of S4 class 'FMM'
fit <- fitFMM(fmm2.data$y,timePoints = fmm2.data$t,nback = 2,
              lengthAlphaGrid = 24,lengthOmegaGrid = 10)
## estimate peakes and trough times and signals
FMMPeaks(fit, timePointsIn2pi = TRUE) # times in the [0,2*pi] interval

}
