\name{hermite.h.recurrences}
\alias{hermite.h.recurrences}
\title{ Recurrence relations for Hermite polynomials }
\description{
  This function returns a data frame with $n$+1 rows and four named columns containing the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations for the order $k$ Hermite polynomial, \eqn{H_k \left( x \right)}, and for orders $k$ = 0, 1, ... , $n$.
}
\usage{
hermite.h.recurrences(n, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{normalized}{ boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\value{
  A data frame with the recurrence relation parameters.
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{hermite.h.inner.products}},
}
\examples{
normalized.r <- hermite.h.recurrences( 10, normalized=TRUE )
unnormalized.r <- hermite.h.recurrences( 10, normalized=FALSE )
}
\keyword{ math }
