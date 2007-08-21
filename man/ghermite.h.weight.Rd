\name{ghermite.h.weight}
\alias{ghermite.h.weight}
\title{ Weight function for the generalized Hermite polynomial }
\description{
  This function returns the value of the weight function for the order $k$ generalized Hermite polynomial, \eqn{
H_k^{\left( \mu  \right)} \left( x \right) }.
}
\usage{
ghermite.h.weight(x, mu)
}
\arguments{
  \item{x}{ a numeric vector function argument }
  \item{mu}{ polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -\infty,\infty \right) }.
  The parameter \eqn{\mu} must be greater than -0.5.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- ghermite.h.weight( x, 1 )
}
\keyword{ math }
