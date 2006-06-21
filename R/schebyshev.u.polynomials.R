schebyshev.u.polynomials <- function( n, normalized=FALSE )
{
###
###	This function returns a list with n+1 elements
###	containing the order k shifted Chebyshev polynomials of the second kind, Uk(x),
###	for orders k=0,1,...n
###
###	Parameters
###	n = integer highest polynomial order
###	normalized = boolean value.  if true, the polynomials are normalized
###
	require( polynom )
	recurrences <- schebyshev.u.recurrences( n, normalized )
	polynomials <- orthogonal.polynomials( recurrences )
	return( polynomials )
}
