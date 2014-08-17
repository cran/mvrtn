\name{rtn}
\alias{rtn}
\title{random variates from truncated normal distribution}

\description{
Random variates from the left or right truncated normal
distribution are generated. An efficient algorithm based
on the inverse CDF method is used. This algorithm can be
used to verify empirically the correcntess of the mean
and variance computation done by the other function
in this package, mvtn. 
}
\usage{
rtn(n, zmu = 0, zsig = 1, c = 1.96, side = c("left", "right"))
}
\arguments{
  \item{n}{
number of random variables 
}
  \item{zmu}{
mean parameter for normal distribution
}
  \item{zsig}{
standard deviation parameter
}
  \item{c}{
truncation point
}
  \item{side}{
either left or right
}
}
\details{
The inverse cdf method is used in the case of the right truncated normal 
distribution. In the left truncated case, the negative of the corresponding
right truncated distribution is used.
}
\value{
the random variates, vector of length n
}

\references{
For further details see the vignette.
}
\author{
A. I. McLeod
}


\seealso{
\code{\link{mvtn}}
}
\examples{
#Empirical validation
z<-rtn(10^4, 100, 15, c=120, side="left")
c(mean(z), var(z))
mvtn(100, 15, c=120, side="left")

}
\keyword{ distribution }
