\name{ghermite.h.polynomials}
\alias{ghermite.h.polynomials}
\title{ Create list of generalized Hermite polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ generalized Hermite polynomials, 
  \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, for orders $k$ = 0, 1, ..., $n$.
}
\usage{
ghermite.h.polynomials(n, mu, normalized = FALSE)
}
\arguments{
  \item{n}{ \code{n} integer highest polynomial order }
  \item{mu}{ \code{mu} polynomial parameter }
  \item{normalized}{ \code{normalized} boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\details{
  The parameter \eqn{\mu} must be greater than -0.5.
}
\value{
  A list of $n$+1 polynomial objects
  \item{1 }{order 0 generalized Hermite polynomial}
  \item{2 }{order 1 generalized Hermite polynomial}
  ...
  \item{n+1 }{order $n$ generalized Hermite polynomial}
}
\references{ Alvarez-Nordase, Atakishiyeva, Atakishiyeva (2004) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{ghermite.h.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
normalized.p.list <- ghermite.h.polynomials( 10, 1, normalized=TRUE )
unnormalized.p.list <- ghermite.h.polynomials( 10, 1, normalized=FALSE )
}
\keyword{ math }
