\name{jacobi.p.inner.products}
\alias{jacobi.p.inner.products}
\title{ Inner products of Jacobi polynomials }
\description{
  This function returns a vector with $n$+1 elements containing the inner product of an order $k$ Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}, with itself (i.e. the norm squared) for orders $k$ = 0, 1, ... , $n$.
}
\usage{
jacobi.p.inner.products(n,a,b)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{a}{ first polynomial parameter }
  \item{b}{ first polynomial parameter }
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
h <- jacobi.p.inner.products( 10, 2, 2 )
}
\keyword{ math }
