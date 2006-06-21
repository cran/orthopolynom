jacobi.p.polynomials <- function( n, a, b, normalized=FALSE )
{
###
###	This function returns a list with n+1 elements
###	containing the order k Jacobi polynomials Pk(a,b,x)
###	for orders k=0,1,...n
###
###	Parameters
###	n = integer highest polynomial order
###	a = first polynomial parameter
###	b = second polynomial parameter
###	normalized = boolean value.  if true, the polynomials are normalized
###
	require( polynom )
	recurrences <- jacobi.p.recurrences( n, a, b, normalized )
	polynomials <- orthogonal.polynomials( recurrences )
	return( polynomials )
}
