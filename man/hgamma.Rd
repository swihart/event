\name{hgamma}
\title{Hazard Function for a Gamma Process}
\alias{hgamma}
\usage{
hgamma(y, shape, rate=1, scale=1/rate)
}
\arguments{
\item{y}{Vector of times.}
\item{shape}{Shape parameter.}
\item{rate}{Rate parameter.}
\item{scale}{Scale parameter.}
}
\value{
\code{hgamma} returns the log hazard function for a gamma process with
the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
\code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
\code{\link[event]{hcauchy}}, \code{\link[rmutil]{hggamma}},
\code{\link[rmutil]{hhjorth}}, \code{\link[rmutil]{hinvgauss}},
\code{\link[rmutil]{hlaplace}}, \code{\link[event]{hlnorm}},
\code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
\code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
\code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
}
\examples{
hgamma(1:10, 3, 2)
}
\keyword{distribution}
