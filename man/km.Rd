\name{km}
\title{Kaplan-Meier Survivor Curves}
\usage{
km(times, censor, group=1, freq=1, cdf=FALSE)
plot.km(z, add=FALSE, xlim, ylim=c(0,1), main=NULL, xlab="Time",
	ylab=NULL, lty=NULL, ...)
plot.intensity.km(z, add=FALSE, xlab="Time", ylab="Intensity", type="l",
	lty=NULL, ...)
plot.dist.km(z)
}
\alias{km}
\alias{plot.km}
\alias{plot.intensity.km}
\alias{plot.dist}
\alias{plot.dist.km}
\description{
\code{km} calculates the Kaplan-Meier estimates for survival.

To plot the survivor curve, use \code{plot()}; for the empirical
intensity curve, use \code{plot.intensity()}; for diagnostic curves to
choose a distribution to which the data might belong, use \code{plot.dist()}.
}
\arguments{
\item{times}{Vector of times to events or a list of vectors of such
times for different individuals.}
\item{censor}{Vector of censoring indicators corresponding to the
vector of times or to the last time in each vector of a list.}
\item{group}{Vector indicating to which group each individual belongs.}
\item{freq}{Vector of frequencies for grouped data.}
\item{cdf}{If TRUE, calculate the cdf instead of the survivor curve.}
\item{z}{An object produced by \code{km}.}
}
\value{
A matrix with class, \code{km}, containing the Kaplan-Meier estimates
is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[event]{plot.intensity}}, \code{\link[event]{plot.survivor}}
}
\examples{
surv <- rgamma(40,2,scale=5)
cens <- rbinom(40,1,0.9)
treat <- gl(2,20)
plot(km(surv, cens, group=treat), main="",xlab="Months",
	ylab="Probability of deterioration")
plot.dist(km(surv, cens, group=treat))
plot.intensity(km(surv, cens, group=treat),ylab="Risk of deterioration")
}
\keyword{hplot}
