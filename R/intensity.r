#
#  event : A Library of Special Functions for Event Histories
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#	hboxcox(y,m,s,f)
#	hburr(y,m,s,f)
#	hcauchy(y,m,s)
#	hexp(y,m)
#	hgextval(y,s,m,f)
#	hgamma(y,s,m)
#	hggamma(y,s,m,f)
#	hhjorth(y,m,s,f)
#	hinvgauss(y,m,s)
#	hlaplace(y,m,s)
#	hlnorm(y,m,s)
#	hlogis(y,m,s)
#	hglogis(y,m,s,f)
#	hnorm(y,m,s)
#	hpareto(y,m,s)
#	hstudent(y,m,s,f)
#	hweibull(y,s,m)
#	hgweibull(y,s,m,f)
#	hskewlaplace(y,s,m,f)
#
#  DESCRIPTION
#
#    Functions for various log hazards or intensities

# for log distributions, subtract log(y) from the intensity function

### Box-Cox intensity
###
# f=1 gives truncated normal
hboxcox <- function(y,m,s,f) {
	y1 <- y^f/f
	-(y1-m)^2/s/2+(f-1)*log(y)-log(2*pi*s)/2-log(1-pnorm(y1,m,sqrt(s))+(f<0)*(1-pnorm(0,m,sqrt(s))))}

### Burr intensity
###
hburr <- function(y,m,s,f) {
	y1 <- y/m
	y2 <- y1^s
	log(f*s/m)+(s-1)*log(y1)-log(1+y2)}

### Cauchy intensity
###








#' Hazard Function for a Cauchy Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param m Location parameter.
#' @param s Dispersion parameter.
#' @return \code{hcauchy} returns the log hazard function for a Cauchy process
#' with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hgamma}}, \code{\link[rmutil]{hggamma}},
#' \code{\link[rmutil]{hhjorth}}, \code{\link[rmutil]{hinvgauss}},
#' \code{\link[rmutil]{hlaplace}}, \code{\link[event]{hlnorm}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hcauchy(1:10, 3, 2)
#' 
#' @export hcauchy
hcauchy <- function(y,m,s) log(dcauchy(y,m,s))-log(1-pcauchy(y,m,s))

### exponential intensity
###








#' Hazard Function for a Poisson Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param rate Vector of rates.
#' @return \code{hexp} returns the log hazard function for a Poisson process
#' with the given parameter value.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hgamma}}, \code{\link[rmutil]{hggamma}},
#' \code{\link[rmutil]{hhjorth}}, \code{\link[rmutil]{hinvgauss}},
#' \code{\link[rmutil]{hlaplace}}, \code{\link[event]{hlnorm}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hexp(1:10, 3)
#' 
#' @export hexp
hexp <- function(y,rate)
	if(length(rate)==1)rep(log(rate),length(y)) else log(rate)

### generalized extreme value intensity
###
# f=1 gives truncated extreme value
hgextval <- function(y,s,m,f) {
	y1 <- y^f/f
	ey <-exp(y1)
	log(s)+s*(y1-log(m))-(ey/m)^s+(f-1)*log(y)-log(1-pweibull(ey,s,m)-(f<0)*exp(-m^-s))}

### gamma intensity
###








#' Hazard Function for a Gamma Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @param scale Scale parameter.
#' @return \code{hgamma} returns the log hazard function for a gamma process
#' with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hcauchy}}, \code{\link[rmutil]{hggamma}},
#' \code{\link[rmutil]{hhjorth}}, \code{\link[rmutil]{hinvgauss}},
#' \code{\link[rmutil]{hlaplace}}, \code{\link[event]{hlnorm}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hgamma(1:10, 3, 2)
#' 
#' @export hgamma
hgamma <- function(y,shape,rate=1,scale=1/rate)
	dgamma(y,shape,scale=scale,log=TRUE)-
		pgamma(y,shape,scale=scale,lower=FALSE,log=TRUE)

### generalized gamma intensity
###
hggamma <- function(y,s,m,f) {
	t <- m/s
	u <- t^f
	y1 <- y^f
	v <- s*f
	-v*log(t)-y1/u+log(f)+(v-1)*log(y)-lgamma(s)-
		log(1-pgamma(y1,s,scale=u))}

### Hjorth intensity
###
hhjorth <- function(y,m,s,f) log(y/m^2+f/(1+s*y))

### inverse Gaussian intensity
###
hinvgauss <- function(y,m,s) {
	t <- y/m
	v <- sqrt(y*s)
	-((t-1)^2/(y*s)+log(2*s*pi*y^3))/2-log(1-pnorm((t-1)/v)
		-exp(2/(m*s))*pnorm(-(t+1)/v))}
### Laplace intensity
###
hlaplace <- function(y,m=0,s=1){
	plp <- function(u){
		t <- exp(-abs(u))/2
		ifelse(u<0,t,1-t)}
	-abs(y-m)/s-log(2*s)-log(1-plp((y-m)/s))}

### log normal intensity
###








#' Hazard Function for a Log Normal Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param m Mean parameter.
#' @param s Variance parameter.
#' @return \code{hlnorm} returns the log hazard function for a log normal
#' process with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hcauchy}}, \code{\link[event]{hgamma}},
#' \code{\link[rmutil]{hggamma}}, \code{\link[rmutil]{hhjorth}},
#' \code{\link[rmutil]{hinvgauss}}, \code{\link[rmutil]{hlaplace}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hnorm}}, \code{\link[event]{hstudent}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hlnorm(1:10, 3, 2)
#' 
#' @export hlnorm
hlnorm <- function(y,m,s) dlnorm(y,m,s,TRUE)-plnorm(y,m,s,FALSE,TRUE)

### logistic intensity
###








#' Hazard Function for a Logistic Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param m Location parameter.
#' @param s Scale parameter.
#' @return \code{hlogis} returns the log hazard function for a logistic process
#' with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link{hboxcox}}, \code{\link{hburr}}, \code{\link{hexp}},
#' \code{\link{hgextval}}, \code{\link{hcauchy}}, \code{\link{hgamma}},
#' \code{\link{hggamma}}, \code{\link{hhjorth}}, \code{\link{hinvgauss}},
#' \code{\link{hlaplace}}, \code{\link{hlnorm}}, \code{\link{hglogis}},
#' \code{\link{hnorm}}, \code{\link{hstudent}}, \code{\link{hweibull}},
#' \code{\link{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hlogis(1:10, 3, 2)
#' 
#' @export hlogis
hlogis <- function(y,m,s) dlogis(y,m,s,TRUE)-plogis(y,m,s,FALSE,TRUE)

### generalized logistic intensity
###
# f=1 gives hlogis
hglogis <- function(y,m,s,f) {
	y1 <- (y-m)/s
	ey <- exp(-y1)
	-log(s/f)-y1-(f+1)*log(1+ey)-log(1-(1+ey)^-f)}

### normal intensity
###








#' Hazard Function for a Normal Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param m Mean parameter.
#' @param s Variance parameter.
#' @return \code{hnorm} returns the log hazard function for a normal process
#' with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hcauchy}}, \code{\link[event]{hgamma}},
#' \code{\link[rmutil]{hggamma}}, \code{\link[rmutil]{hhjorth}},
#' \code{\link[rmutil]{hinvgauss}}, \code{\link[rmutil]{hlaplace}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hlnorm}}, \code{\link[event]{hstudent}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hnorm(1:10, 3, 2)
#' 
#' @export hnorm
hnorm <- function(y,m,s) dnorm(y,m,s,TRUE)-pnorm(y,m,s,FALSE,TRUE)

### Pareto intensity
###
hpareto <- function(y,m,s) (s+1)/(m*s+y)

### Student t intensity
###








#' Hazard Function for a Student t Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param m Location parameter.
#' @param s Scale parameter.
#' @param f Degrees of freedom.
#' @return \code{hstudent} returns the log hazard function for a Student t
#' process with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hcauchy}}, \code{\link[event]{hgamma}},
#' \code{\link[rmutil]{hggamma}}, \code{\link[rmutil]{hhjorth}},
#' \code{\link[rmutil]{hinvgauss}}, \code{\link[rmutil]{hlaplace}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hnorm}}, \code{\link[event]{hnorm}},
#' \code{\link[event]{hweibull}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hstudent(1:10, 3, 2, 5)
#' 
#' @export hstudent
hstudent <- function(y,m,s,f){
	pst <- function(u,f){
		t <- 0.5*pbeta(f/(f+u^2),f/2,0.5)
		ifelse(u<0,t,1-t)}
	t <- (f+1)/2
	u <- (y-m)/s
	lgamma(t)-lgamma(f/2)-log(f)/2-(t)*log(1+u^2/f)
		-log(pi)/2-log(1-pst(u,f))}
### Weibull intensity
###








#' Hazard Function for a Weibull Process
#' 
#' 
#' 
#' @param y Vector of times.
#' @param s Shape parameter.
#' @param m Scale parameter.
#' @return \code{hweibull} returns the log hazard function for a Weibull
#' process with the given parameter values.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{hboxcox}}, \code{\link[rmutil]{hburr}},
#' \code{\link[event]{hexp}}, \code{\link[rmutil]{hgextval}},
#' \code{\link[event]{hcauchy}}, \code{\link[event]{hgamma}},
#' \code{\link[rmutil]{hggamma}}, \code{\link[rmutil]{hhjorth}},
#' \code{\link[rmutil]{hinvgauss}}, \code{\link[rmutil]{hlaplace}},
#' \code{\link[event]{hlogis}}, \code{\link[rmutil]{hglogis}},
#' \code{\link[event]{hlnorm}}, \code{\link[event]{hnorm}},
#' \code{\link[event]{hstudent}}, \code{\link[rmutil]{hgweibull}}.
#' @keywords distribution
#' @examples
#' 
#' hweibull(1:10, 1.5, 2)
#' 
#' @export hweibull
hweibull <- function(y,s,m) log(s)+(s-1)*log(y)-s*log(m)

### generalized Weibull intensity
###
# Mudholkar, Srivastava, & Freimer (1995) Technometrics 37: 436-445
hgweibull <- function(y,s,m,f) {
	y1 <- y/m
	y2 <- y1^s
	y3 <- exp(-y2)
	log(s*f/m)+(s-1)*log(y1)+(f-1)*log(1-y3)-y2-log(1-(1-y3)^f)}

### skew Laplace intensity
###
hskewlaplace <- function(y,m=0,s=1,f=1){
	plp <- function(u)
		ifelse(u>0,1-exp(-f*abs(u))/(1+f^2),f^2*exp(-abs(u)/f)/(1+f^2))
	log(f)+ifelse(y>m,-f*(y-m),(y-m)/f)/s-log((1+f^2)*s)-
		log(1-plp((y-m)/s))}
