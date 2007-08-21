\name{ghermite.h.inner.products}
\alias{ghermite.h.inner.products}
\title{ Inner products of generalized Hermite polynomials }
\description{
  This function returns a vector with $n$+1 elements containing the inner product of an order $k$ generalized Hermite polynomial,
   \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, with itself (i.e. the norm squared) for orders $k$ = 0, 1, ... , $n$.
}
\usage{
ghermite.h.inner.products(n, mu)
}
\arguments{
  \item{n}{ \code{n} integer highest polynomial order }
  \item{mu}{ \code{mu} polynomial parameter }
}
\details{
  The parameter \eqn{\mu} must be greater than -0.5.
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
h <- ghermite.h.inner.products( 10, 1 )
}
\keyword{ math }
