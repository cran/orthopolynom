\name{legendre.weight}
\alias{legendre.weight}
\title{ Weight function for the Legendre polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Legendre polynomial, \eqn{P_k \left( x \right)}.
}
\usage{
legendre.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -1,1 \right) }.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- legendre.weight( x )
}
\keyword{ math }
