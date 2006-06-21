\name{polynomial.coefficients}
\alias{polynomial.coefficients}
\title{ Create list of polynomial coefficient vectors }
\description{
  This function returns a list with $n$+1 elements containing the vector of coefficients of the order $k$ polynomials for orders $k$ = 0, 1,..., $n$.  Each element in the list is a vector.
}
\usage{
polynomial.coefficients(polynomials)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{polynomials}{ list of polynomial objects }
}
\value{
  A list of $n$+1 polynomial objects where each element is a vector of coefficients.
  \item{1 }{Coefficient(s) of the order 0 polynomial}
  \item{2 }{Coefficient(s) of the order 1 polynomial}
  ...
  \item{n+1}{Coefficient(s) of the order $n$ polynomial}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link[polynom]{predict.polynomial}}
}
\examples{
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
p.coef <- polynomial.coefficients( p.list )
}
\keyword{ math }
