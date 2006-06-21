\name{polynomial.values}
\alias{polynomial.values}
\title{ Create vector of polynomial values }
\description{
  This function returns a list with $n$+1 elements containing the values of the order $k$ polynomials for orders $k$ = 0, 1,..., $n$ and for the given argument $x$.
}
\usage{
polynomial.values( polynomials, x )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{polynomials}{ list of polynomial objects }
  \item{x}{ the argument which can be any numeric object }
}
\value{
  A list of $n$+1 polynomial objects where each element is the value of the polynomial.
  \item{1 }{Value(s) for the order 0 polynomial}
  \item{2 }{Value(s) for the order 1 polynomial}
  ...
  \item{n+1}{Value(s) for the order $n$ polynomial}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link[polynom]{predict.polynomial}}
}
\examples{
p.list <- chebyshev.t.polynomials( 10, normalized=FALSE )
x <- seq( -2, 2, .01 )
y <- polynomial.values( p.list, x )
}
\keyword{ math }
