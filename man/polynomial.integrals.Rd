\name{polynomial.integrals}
\alias{polynomial.integrals}
\title{ Create list of polynomial integrals }
\description{
  This function returns a list with $n$+1 elements containing polynomial objects which are the indefinite integrals of the order $k$ polynomials for orders $k$ = 0, ... , $n$.
}
\usage{
polynomial.integrals(polynomials)
}
\arguments{
  \item{polynomials}{ list of polynomial objects }
}
\details{
  The polynomial objects in the argument \code{polynomials} are as follows
  \itemize{
  \item{1}{order 0 polynomial}
  \item{2}{order 1 polynomial}
  ...
  \item{n+1}{order n polynomial}
  }
}
\value{
  List of polynomial objects
  \item{1 }{integral of \code{polynomials[[1]]}}
  \item{2 }{integral of \code{polynomials[[2]]}}
  ...
  \item{n+1 }{integral of \code{polynomials[[n+1]]}}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{ 
  \code{\link[polynom]{integral.polynomial}}
}
\examples{
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
p.int <- polynomial.integrals( p.list )
}
\keyword{ math }
