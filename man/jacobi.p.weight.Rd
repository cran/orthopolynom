\name{jacobi.p.weight}
\alias{jacobi.p.weight}
\title{ Weight function for the Jacobi polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}.
}
\usage{
jacobi.p.weight(x,alpha,beta)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ the first polynomial parameter }
  \item{beta}{ the second polynomial parameter }
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
y <- jacobi.p.weight( x, 2, 2 )
}
\keyword{ math }
