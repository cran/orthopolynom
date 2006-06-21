\name{laguerre.inner.products}
\alias{laguerre.inner.products}
\title{ Inner products of Laguerre polynomials }
\description{
  This function returns a vector with $n$+1 elements containing the inner product of an order $k$ Laguerre polynomial, \eqn{L_n \left( x \right)}, with itself (i.e. the norm squared) for orders $k$ = 0, 1, ... , $n$.
}
\usage{
laguerre.inner.products(n)
}
\arguments{
  \item{n}{ integer highest polynomial order }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  A vector with $n$+1 elements
  \item{1 }{inner product of order 0 orthogonal polynomial}
  \item{2 }{inner product of order 1 orthogonal polynomial}
  ...
  \item{n+1 }{inner product of order $n$+1 orthogonal polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
h <- laguerre.inner.products( 10 )
}
\keyword{ math }
