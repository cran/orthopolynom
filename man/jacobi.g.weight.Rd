\name{jacobi.g.weight}
\alias{jacobi.g.weight}
\title{ Weight function for the Jacobi polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Jacobi polynomial, \eqn{G_k \left( {p,q,x} \right)}.
}
\usage{
jacobi.g.weight(x,p,q)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{p}{ the first polynomial parameter }
  \item{q}{ the second polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,1 \right) }.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -2, 2, .01 )
y <- jacobi.g.weight( x, 3, 2 )
}
\keyword{ math }
