\name{chebyshev.c.weight}
\alias{chebyshev.c.weight}
\title{ Weight function for the Chebyshev polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Chebyshev polynomial of the first kind, \eqn{C_k \left( x \right) }.
}
\usage{
chebyshev.c.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -2,2 \right) }.
}
\value{
  The value of the weight function
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
x <- seq( -3, 3, .01 )
y <- chebyshev.c.weight( x )
}
\keyword{ math }
