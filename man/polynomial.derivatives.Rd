\name{polynomial.derivatives}
\alias{polynomial.derivatives}
\title{ Create list of polynomial derivatives }
\description{
  This function returns a list with $n$+1 elements containing polynomial objects which are the derivatives of the order $k$ polynomials for orders $k$ = 0, ... , $n$.
}
\usage{
polynomial.derivatives(polynomials)
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
  \item{1 }{derivative of \code{polynomials[[1]]}}
  \item{2 }{derivative of \code{polynomials[[2]]}}
  ...
  \item{n+1 }{derivative of \code{polynomials[[n+1]]}}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{ 
  \code{\link[polynom]{deriv.polynomial}}
}
\examples{
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
p.deriv <- polynomial.derivatives( p.list )
}
\keyword{ math }
