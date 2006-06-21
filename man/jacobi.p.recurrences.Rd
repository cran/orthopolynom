\name{jacobi.p.recurrences}
\alias{jacobi.p.recurrences}
\title{ Recurrence relations for Jacobi polynomials }
\description{
  This function returns a data frame with $n$+1 rows and four named columns containing the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations for the order $k$ Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}, and for orders $k$ = 0, 1, ... , $n$.
}
\usage{
jacobi.p.recurrences(n, a, b, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{a}{ first polynomial parameter }
  \item{b}{ second polynomial parameter }
  \item{normalized}{ boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\value{
  A data frame with the recurrence relation parameters.
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{jacobi.p.inner.products}},
\code{\link{pochhammer}}
}
\examples{
normalized.r <- jacobi.p.recurrences( 10, 2, 2, normalized=TRUE )
unnormalized.r <- jacobi.p.recurrences( 10, 2, 2, normalized=FALSE )
}
\keyword{ math }
