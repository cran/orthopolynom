\name{schebyshev.u.weight}
\alias{schebyshev.u.weight}
\title{ Weight function for the shifted Chebyshev polynomial }
\description{
  This function returns the value of the weight function for the order $k$ shifted Chebyshev polynomial of the second kind, \eqn{U_k^* \left( x \right) }.
}
\usage{
schebyshev.u.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,1 \right) }.
}
\value{
  The value of the weight function.
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- schebyshev.u.weight( x )
}
\keyword{ math }