% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betasim_funs.R
\name{rtrpois}
\alias{dtrpois}
\alias{rtrpois}
\title{Truncated Poisson distribution}
\usage{
rtrpois(n, lambda)

dtrpois(x, lambda, log = FALSE)
}
\arguments{
\item{n}{number of random deviations}

\item{lambda}{mean of \emph{truncated} distribution}

\item{x}{observed number of counts}

\item{log}{return log-density?}
}
\description{
Truncated Poisson distribution
}
\examples{
set.seed(101)
pp <- prop.table(table(rtrpois(100000,lambda=4.5)))
par(las=1,bty="l")
plot(pp,type="h")
points(1:15,dtrpois(1:15,4.5),col="red",pch=16,cex=2)
}

