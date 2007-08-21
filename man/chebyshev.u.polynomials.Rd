\name{chebyshev.u.polynomials}
\alias{chebyshev.u.polynomials}
\title{ Create list of Chebyshev polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ Chebyshev polynomials of the second kind, \eqn{U_k \left( x\right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
chebyshev.u.polynomials(n, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function produces a data frame with the recurrence relation parameters for the orthogonal polynomials.  It then uses the function \code{orthogonal.polynomials} to construct the list of polynomial objects from the recurrence relations.
}
\value{
  A list of $n$+1 polynomial objects
  \item{1 }{order 0 Chebyshev polynomial}
  \item{2 }{order 1 Chebyshev polynomial}
  ...
  \item{n+1 }{order $n$ Chebyshev polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{chebyshev.u.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- chebyshev.u.polynomials( 10, normalized=TRUE )
unnormalized.p.list <- chebyshev.u.polynomials( 10, normalized=FALSE )
}
\keyword{ math }
