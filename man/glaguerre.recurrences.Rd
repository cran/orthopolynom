\name{glaguerre.recurrences}
\alias{glaguerre.recurrences}
\title{ Recurrence relations for generalized Laguerre polynomials }
\description{
  This function returns a data frame with $n$+1 rows and four named columns containing the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations for the order $k$ generalized Laguerre polynomial, \eqn{L_n^{\left( \alpha  \right)} \left( x \right)}, and for orders $k$ = 0, 1, ... , $n$.
}
\usage{
glaguerre.recurrences(n, alpha, normalized=FALSE)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{alpha}{ polynomial parameter }
  \item{normalized}{ boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\value{
  A data frame with the recurrence relation parameters.
}
\references{ Abramowitz and Stegun (1968) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{glaguerre.inner.products}},
}
\examples{
normalized.r <- glaguerre.recurrences( 10, 1, normalized=TRUE )
unnormalized.r <- glaguerre.recurrences( 10, 1, normalized=FALSE )
}
\keyword{ math }
