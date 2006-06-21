\name{glaguerre.weight}
\alias{glaguerre.weight}
\title{ Weight function for the generalized Laguerre polynomial }
\description{
  This function returns the value of the weight function for the order $k$ generalized Laguerre polynomial, \eqn{L_n^{\left( \alpha  \right)} \left( x \right)}.
}
\usage{
glaguerre.weight(x,alpha)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ polynomial parameter }
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
y <- glaguerre.weight( x, 1 )
}
\keyword{ math }
