\name{chebyshev.s.inner.products}
\alias{chebyshev.s.inner.products}
\title{ Inner products of Chebyshev polynomials }
\description{
  This function returns a vector with $n$+1 elements containing the inner product of an order $k$ Chebyshev polynomial of the second kind, \eqn{S_k \left( x\right)}, with itself (i.e. the norm squared) for orders $k$ = 0, 1, ... , $n$.
}
\usage{
chebyshev.s.inner.products(n)
}
%- maybe also 'usage' for other objects documented here.
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
h <- chebyshev.s.inner.products( 10 )
}
\keyword{ math }
