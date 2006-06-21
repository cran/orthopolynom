\name{chebyshev.u.recurrences}
\alias{chebyshev.u.recurrences}
\title{ Recurrence relations for Chebyshev polynomials }
\description{
  This function returns a data frame with $n$+1 rows and four named columns containing the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations for the order $k$ Chebyshev polynomial of the second kind, \eqn{U_k \left( x \right)}, and for orders $k$ = 0, 1, ... , $n$.
}
\usage{
chebyshev.u.recurrences(n, normalized=FALSE)
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
\code{\link{chebyshev.u.inner.products}},
}
\examples{
normalized.r <- chebyshev.u.recurrences( 10, normalized=TRUE )
unnormalized.r <- chebyshev.u.recurrences( 10, normalized=FALSE )
}
\keyword{ math }