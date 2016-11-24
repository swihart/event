\name{hlnorm}
\title{Hazard Function for a Log Normal Process}
\alias{hlnorm}
\usage{
hlnorm(y, m, s)
}
\arguments{
\item{y}{Vector of times.}
\item{m}{Mean parameter.}
\item{s}{Variance parameter.}
}
\value{
\code{hlnorm} returns the log hazard function for a log normal
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
\code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
\code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
}
\examples{
hlnorm(1:10, 3, 2)
}
\keyword{distribution}
