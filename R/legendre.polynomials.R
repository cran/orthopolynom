legendre.polynomials <- function( n, normalized=FALSE )
{
###
###   This function returns a list with n+1 elements
###   containing the order k Legendre polynomials Pk(x),
###   for orders k=0,1,...,n
###
###   Parameters
###   n = integer highest polynomial order
###   normalized = boolean value.  if true, the polynomials are normalized
###
   require( polynom )
   recurrences <- legendre.recurrences( n, normalized )
   polynomials <- orthogonal.polynomials( recurrences )
   return( polynomials )
}
