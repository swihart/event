\name{hstudent}
\title{Hazard Function for a Student t Process}
\alias{hstudent}
\usage{
hstudent(y, m, s, f)
}
\arguments{
\item{y}{Vector of times.}
\item{m}{Location parameter.}
\item{s}{Scale  parameter.}
\item{f}{Degrees of freedom.}
}
\value{
\code{hstudent} returns the log hazard function for a Student t
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
\code{\link[event]{hnorm}}, \code{\link[event]{hnorm}},
\code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
}
\examples{
hstudent(1:10, 3, 2, 5)
}
\keyword{distribution}
