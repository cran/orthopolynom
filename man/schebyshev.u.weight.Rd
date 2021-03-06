\name{schebyshev.u.weight}
\alias{schebyshev.u.weight}
\title{ Weight function for the shifted Chebyshev polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  shifted Chebyshev polynomial of the second kind, \eqn{U_k^* \left( x \right) }.
}
\usage{
schebyshev.u.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,1 \right) }.  The formula
  used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \sqrt {x - x^2 } }
}
\value{
  The value of the weight function.
}
\references{
Abramowitz, M. and I. A. Stegun, 1968. \emph{Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., New York.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
###
### compute the shifted U Chebyshev weight function for argument values
### between 0 and 1
###
x <- seq( 0, 1, .01 )
y <- schebyshev.u.weight( x )
plot( x, y )
}
\keyword{ math }
