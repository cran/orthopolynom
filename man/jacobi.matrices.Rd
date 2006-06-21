\name{jacobi.matrices}
\alias{jacobi.matrices}
\title{ Create list of Jacobi matrices from monic recurrence parameters }
\description{
  Return a list of $n$ real symmetric, tri-diagonal matrices which are the principal minors of the \eqn{n \times n} Jacobi matrix derived from the monic recurrence parameters, $a$ and $b$, for orthogonal polynomials.
}
\usage{
jacobi.matrices(r)
}
\arguments{
  \item{r}{ a data frame containing the parameters $a$ and $b$ }
}
\value{
  A list of symmetric, tri-diagnonal matrices
  \item{1 }{a \eqn{1 \times 1} matrix}
  \item{2 }{a \eqn{2 \times 2} matrix}
  ...
  \item{n }{an \eqn{n \times n} matrix}
}
\references{ Press, et. al. (1992) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
r <- chebyshev.t.recurrences( 5 )
m.r <- monic.polynomial.recurrences( r )
j.m <- jacobi.matrices( m.r )
}
\keyword{ math }
