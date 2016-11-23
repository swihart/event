\name{hweibull}
\title{Hazard Function for a Weibull Process}
\alias{hweibull}
\usage{
hweibull(y, s, m)
}
\arguments{
\item{y}{Vector of times.}
\item{s}{Shape parameter.}
\item{m}{Scale parameter.}
}
\value{
\code{hweibull} returns the log hazard function for a Weibull
process with the given parameter values.
}
\author{J.K. Lindsey}
\seealso{
\code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
\code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
\code{\link[event]{hcauchy}}, \code{\link[event]{hgamma}},
\code{\link[rmutil]{hggamma}}, \code{\link[rmutil]{hhjorth}},
\code{\link[rmutil]{hinvgauss}}, \code{\link[rmutil]{hlaplace}},
\code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
\code{\link[event]{hlnorm}}, \code{\link[event]{hnorm}},
\code{\link[event]{hstudent}}, \code{\link[rmutil]{hgweibull}}.
}
\examples{
hweibull(1:10, 1.5, 2)
}
\keyword{distribution}
