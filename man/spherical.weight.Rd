\name{spherical.weight}
\alias{spherical.weight}
\title{ Weight function for the spherical polynomial }
\description{
  This function returns the value of the weight function for the order $k$ spherical polynomial, \eqn{P_k \left( x \right)}.
}
\usage{
spherical.weight(x)
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
y <- spherical.weight( x )
}
\keyword{ math }
