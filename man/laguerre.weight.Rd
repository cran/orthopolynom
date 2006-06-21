\name{laguerre.weight}
\alias{laguerre.weight}
\title{ Weight function for the Laguerre polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Laguerre polynomial, \eqn{L_n \left( x \right)}.
}
\usage{
laguerre.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,\infty \right) }.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- laguerre.weight( x )
}
\keyword{ math }
