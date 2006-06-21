\name{polynomial.orders}
\alias{polynomial.orders}
\title{ Create vector of polynomial orders }
\description{
  This function returns a vector with $n$ elements containing the orders of the polynomials
}
\usage{
polynomial.orders(polynomials)
}
\arguments{
  \item{polynomials}{ list of $n$ polynomial objects }
}
\value{
  A vector of $n$ values
  \item{1 }{Order of polynomials[[1]]}
  \item{2 }{Order of polynomials[[2]]}
  ...
  \item{n}{Order of polynomials[[n]]}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
p.order <- polynomial.orders( p.list )
}
\keyword{ math }
