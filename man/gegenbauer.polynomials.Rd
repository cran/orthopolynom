\name{gegenbauer.polynomials}
\alias{gegenbauer.polynomials}
\title{ Create list of Gegenbauer polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ Gegenbauer polynomials, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
gegenbauer.polynomials(n, alpha, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{alpha}{ polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function produces a data frame with the recurrence relation parameters for the orthogonal polynomials.  It then uses the function \code{orthogonal.polynomials} to construct the list of polynomial objects from the recurrence relations.
}
\value{
  A list of $n$+1 polynomial objects
  \item{1 }{order 0 Gegenbauer polynomial}
  \item{2 }{order 1 Gegenbauer polynomial}
  ...
  \item{n+1 }{order $n$ Chebyshev polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{gegenbauer.recurrences}},
\code{\link{orthogonal.polynomials}}
}
\examples{
normalized.p.list <- gegenbauer.polynomials( 10, 1, normalized=TRUE )
unnormalized.p.list <- gegenbauer.polynomials( 10, 1, normalized=FALSE )
}
\keyword{ math }
