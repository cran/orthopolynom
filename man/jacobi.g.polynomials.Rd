\name{jacobi.g.polynomials}
\alias{jacobi.g.polynomials}
\title{ Create list of Jacobi polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ Jacobi polynomials, \eqn{G_k \left( {p,q,x} \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
jacobi.g.polynomials(n, p, q, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{p}{ first polynomial parameter }
  \item{q}{ second polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function produces a data frame with the recurrence relation parameters for the orthogonal polynomials.  It then uses the function \code{orthogonal.polynomials} to construct the list of polynomial objects from the recurrence relations.
}
\value{
  A list of $n$+1 polynomial objects
  \item{1 }{order 0 Jacobi polynomial}
  \item{2 }{order 1 Jacobi polynomial}
  ...
  \item{n+1 }{order $n$ Chebyshev polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{jacobi.g.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- jacobi.g.polynomials( 10, 3, 2, normalized=TRUE )
unnormalized.p.list <- jacobi.g.polynomials( 10, 3, 2, normalized=FALSE )
}
\keyword{ math }
