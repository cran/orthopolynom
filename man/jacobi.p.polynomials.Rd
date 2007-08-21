\name{jacobi.p.polynomials}
\alias{jacobi.p.polynomials}
\title{ Create list of Jacobi polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ Jacobi polynomials, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
jacobi.p.polynomials(n, a, b, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{a}{ first polynomial parameter }
  \item{b}{ second polynomial parameter }
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
\code{\link{jacobi.p.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- jacobi.p.polynomials( 10, 2, 2, normalized=TRUE )
unnormalized.p.list <- jacobi.p.polynomials( 10, 2, 2, normalized=FALSE )
}
\keyword{ math }
