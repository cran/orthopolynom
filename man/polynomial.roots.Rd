\name{polynomial.roots}
\alias{polynomial.roots}
\title{ Create a list of polynomial roots }
\description{
  This function returns a list with $n$ elements containing the roots of the order $k$ monic orthogonal polynomials for orders $k$ = 1, ... , $n$ using a data frame with the monic polynomial recurrence parameter vectors \eqn{\bf{a}} and \eqn{\bf{b}}
}
\usage{
polynomial.roots(m.r)
}
\arguments{
  \item{m.r}{ monic recurrence data frame with parameters a and b }
}
\details{
  The parameter \code{m.r} is a data frame with $n$+1 rows and two names columns.  The columns which are names \code{a} and \code{b} correspond to the above referenced vectors.
  Function \code{jacobi.matrices} is used to create a list of symmetric, tridiagonal Jacobi matrices from these named columns.  The eigenvalues of the \eqn{k \times k} Jacobi matrix are the roots or zeros of the order $k$ monic orthogonal polynomial.
}
\value{
  A list with n elements each of which is a vector of polynomial roots
  \item{1 }{roots of the order 1 monic polynomial}
  \item{2 }{roots of the order 2 monic polynomial}
  ...
  \item{n }{roots of the order $n$ monic polynomial}
}
\references{ Press, et. al. (1992) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{ 
\code{\link{monic.polynomial.recurrences}}
\code{\link{jacobi.matrices}},
\code{\link[base]{eigen}}
}
\examples{
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
m.r <- monic.polynomial.recurrences( r )
p.roots <- polynomial.roots( m.r )
}
\keyword{ math }
