\name{monic.polynomials}
\alias{monic.polynomials}
\title{ Create list of monic orthogonal polynomials }
\description{
  This function returns a list with $n$+1 elements containing the order $k$ monic polynomials for orders $k$ = 0, 1 ,... ,$n$.
}
\usage{
monic.polynomials(monic.recurrences)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monic.recurrences}{ a data frame containing the parameters a and b }
}
\value{
  A list with $n$+1 polynomial objects
  \item{1 }{order 0 monic orthogonal polynomial}
  \item{2 }{order 1 monic orthogonal polynomial}
  ...
  \item{n+1 }{order $n$ monic orthogonal polynomial}
}
\references{ Press, et al (1992) }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{monic.polynomial.recurrences}}
}
\examples{
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
m.r <- monic.polynomial.recurrences( r )
p.list <- monic.polynomials( m.r )
}
\keyword{ math }
