\name{orthonormal.polynomials}
\alias{orthonormal.polynomials}
\title{ Create orthonormal polynomials }
\description{
  Create list of orthonormal polynomials from  the following recurrence relations for $k$ = 0, 1,...$n$.
  
  \deqn{c_k p_{k+1}\left( x \right) = \left( d_k + e_k x \right) p_k \left( x \right) - f_k p_{k-1} \left( x \right)}
  
  We require that \eqn{p_{-1} \left( x \right) = 0} and \eqn{p_0 \left( x \right) = 1}.  The coefficients are the column vectors \eqn{{\bf{c}}}, \eqn{{\bf{d}}}, \eqn{{\bf{e}}} and \eqn{{\bf{f}}}.
}
\usage{
orthonormal.polynomials(recurrences, p.0)
}
\arguments{
  \item{recurrences}{ a data frame containing the parameters of the orthonormal polynomial recurrence relations }
  \item{p.0}{ a polynomial object for the order 0 orthonormal polynomial }
}
\details{
  The argument is a data frame with $n$+1 rows and four named columns.  The column names are \code{c}, \code{d}, \code{e} and \code{f}.  These columns correspond to the column vectors described above.
}
\value{
  A list of polynomial objects
  \item{1 }{Order 0 orthonormal polynomial}
  \item{2 }{Order 1 orthonormal polynomial}
  ...
  \item{n+1}{Order $n$ orthonormal polynomial}
}
\references{ Abramowitz and Stegun (1968) and Press, et. al. (1992)}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
norm <- sqrt( pi )
p.0 <- polynomial( c( 1 / norm ) )
p.list <- orthonormal.polynomials( r, p.0 )
}
\keyword{ math }
