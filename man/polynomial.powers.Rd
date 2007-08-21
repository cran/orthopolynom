\name{polynomial.powers}
\alias{polynomial.powers}
\title{ Create a list of polynomial linear combinations }
\description{
  This function returns a list with $n$+1 elements containing the vector of 
  linear combinations of the order $k$ polynomials for orders $k$ = 0, 1,..., $n$.
  Each element in the list is a vector.
}
\usage{
polynomial.powers(polynomials)
}
\arguments{
  \item{polynomials}{ ~~Describe \code{polynomials} here~~ }
}
\details{
  The $j$-th component in the list is a vector of linear combinations of
  the order $k$ polynomials for orders $k$ = 0, 1, ..., $j-1$ equal to
  the monomial x raised to the power $j$-1.
}
\value{
  A list of $n$+1 elements where each element is a vector of linear combinations.
  \item{1 }{Linear combination(s) of polynomials up to order 0}
  \item{2 }{Linear combination(s) of polynomials up to order 1}
  ...
  \item{n+1}{Linear combination(s) of polynomials up to order $n$}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
polynomials <- legendre.polynomials( 10 )
alphas <- polynomial.powers( polynomials )
}
\keyword{ math }
