\name{glaguerre.polynomials}
\alias{glaguerre.polynomials}
\title{ Create list of generalized Laguerre polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ generalized Laguerre polynomials, \eqn{L_n^{\left( \alpha  \right)} \left( x \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
glaguerre.polynomials(n, alpha, normalized=FALSE)
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
  \item{1 }{order 0 generalized Laguerre polynomial}
  \item{2 }{order 1 generalized Laguerre polynomial}
  ...
  \item{n+1 }{order $n$ generalized Laguerre polynomial}
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{glaguerre.recurrences}},
\code{\link{orthogonal.polynomials}}
}
\examples{
normalized.p.list <- glaguerre.polynomials( 10, 1, normalized=TRUE )
unnormalized.p.list <- glaguerre.polynomials( 10, 1, normalized=FALSE )
}
\keyword{ math }
