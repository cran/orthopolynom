\name{ultraspherical.inner.products}
\alias{ultraspherical.inner.products}
\title{ Inner products of ultraspherical polynomials }
\description{
  This function returns a vector with $n$+1 elements containing the inner product of an order $k$ ultraspherical polynomial, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}, with itself (i.e. the norm squared) for orders $k$ = 0, 1, ... , $n$.
}
\usage{
ultraspherical.inner.products(n,alpha)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{alpha}{ polynomial parameter }
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
\seealso{
\code{\link{gegenbauer.inner.products}}
}
\examples{
h <- ultraspherical.inner.products( 10, 1 )
}
\keyword{ math }