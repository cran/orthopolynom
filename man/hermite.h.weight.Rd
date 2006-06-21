\name{hermite.h.weight}
\alias{hermite.h.weight}
\title{ Weight function for the Hermite polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Hermite polynomial, \eqn{H_k \left( x \right)}.
}
\usage{
hermite.h.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -\infty,\infty \right) }.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- hermite.h.weight( x )
}
\keyword{ math }
