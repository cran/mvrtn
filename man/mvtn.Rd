\name{mvtn}
\alias{mvtn}

\title{
Mean and Variance of Truncated Normal Distribution
}
\description{
Computes the mean and variance in a univariate censored normal
distribution with parameters zmu and zsig and censor point c.
}
\usage{
mvtn(zmu, zsig, c, side = c("left", "right"))
}

\arguments{
  \item{zmu}{
    Mean Paramter
}
  \item{zsig}{
    Standard Deviation Paramter
}
  \item{c}{
    Censor Point
}
  \item{side}{
  Left or right truncation of normal distribution
}
}

\details{
The derivation of the method used is described in the vignette
that accompanies this package (Mohammad, McLeod, and McLeod, 2014).
See documentation for the accompanying function rtn() for
empirical validation.
}

\value{
a list with two components:
\item{mean}{mean of the distribution}
\item{variance}{variance of the distribution}
}

\references{
N. Mohammad, M. McLeod and I. McLeod, Mean and Variance of the Truncated Normal Distribution.
Vignette for mvtn package.
}

\author{
Matthew McLeod
}

\seealso{
\code{\link{rtn}}
}

\examples{
mvtn(100,15,80,"left")
mvtn(100,15,80,"right")
#simulation check on the answer
\dontrun{
#Simulate truncated left truncated normal and compute
#  empirical means and variances and their sds

}
}

\keyword{ distribution }

