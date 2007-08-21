\name{schebyshev.t.polynomials}
\alias{schebyshev.t.polynomials}
\title{ Create list of shifted Chebyshev polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ shifted Chebyshev polynomials of the first kind, \eqn{T_k^* \left( x\right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
schebyshev.t.polynomials(n, normalized)
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
  \item{1 }{order 0 shifted Chebyshev polynomial}
  \item{2 }{order 1 shifted Chebyshev polynomial}
  ...
  \item{n+1 }{order $n$ shifted Chebyshev polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{schebyshev.u.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- schebyshev.t.polynomials( 10, normalized=TRUE )
unnormalized.p.list <- schebyshev.t.polynomials( 10, normalized=FALSE )
}
\keyword{ math }
