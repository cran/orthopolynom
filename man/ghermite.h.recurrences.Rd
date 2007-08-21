\name{ghermite.h.recurrences}
\alias{ghermite.h.recurrences}
\title{ Recurrence relations for generalized Hermite polynomials }
\description{
  This function returns a data frame with $n$+1 rows and four named columns containing the coefficient vectors
   \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations for the order $k$ generalized Hermite polynomial,
    \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, and for orders $k$ = 0, 1, ... , $n$.
}
\usage{
ghermite.h.recurrences(n, mu, normalized = FALSE)
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
  A data frame with the recurrence relation parameters.
}
\references{ Alvarez-Nordase, Atakishiyeva, Atakishiyeva (2004) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{ghermite.h.inner.products}},
}
\examples{
normalized.r <- ghermite.h.recurrences( 10, 1, normalized=TRUE )
unnormalized.r <- ghermite.h.recurrences( 10, 1, normalized=FALSE )
}
\keyword{ math }
