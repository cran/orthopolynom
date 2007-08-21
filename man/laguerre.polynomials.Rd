\name{laguerre.polynomials}
\alias{laguerre.polynomials}
\title{ Create list of Laguerre polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ Laguerre polynomials, \eqn{L_n \left( x \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
laguerre.polynomials(n, normalized=FALSE)
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
  \item{1 }{order 0 Laguerre polynomial}
  \item{2 }{order 1 Laguerre polynomial}
  ...
  \item{n+1 }{order $n$ Laguerre polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{laguerre.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- laguerre.polynomials( 10, normalized=TRUE )
unnormalized.p.list <- laguerre.polynomials( 10, normalized=FALSE )
}
\keyword{ math }
