\name{gegenbauer.weight}
\alias{gegenbauer.weight}
\title{ Weight function for the Gegenbauer polynomial }
\description{
  This function returns the value of the weight function for the order $k$ Gegenbauer polynomial, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}.
}
\usage{
gegenbauer.weight(x,alpha)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ polynomial parameter }
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
y <- gegenbauer.weight( x, 1 )
}
\keyword{ math }
