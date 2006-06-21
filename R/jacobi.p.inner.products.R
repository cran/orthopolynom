jacobi.p.inner.products <- function( n, a, b )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k Jacobi polynomial
###	Pk(a,b,x) for orders k=0,1,...,n
###
###	Parameters
###	n = integer highest polynomial order
###	a = first parameter
###	b = second parameter
###
	if ( n < 0 )
		stop( "negative highest polynomial order" )
	if ( n != round( n ) )
		stop( "highest polynomial order is not integer" )
	if ( a <= -1 )
		stop( "alpha less than or equal to -1" )
	if ( b <= -1 )
		stop( "beta less than or equal to -1" )
	ab <- a + b	
	coef.num <- 2 ^ ( ab + 1 )
	inner.products <- rep( 1, n + 1 )
	j <- 1
	for ( k in 0:n ) {
		coef.den <- 2 * k + ab + 1
		log.factor <- lgamma( k + a + 1 ) + lgamma( k + b + 1 ) - lfactorial( k ) - lgamma( k + ab + 1 )
		inner.products[j] <- (coef.num / coef.den ) * exp( log.factor )
		j <- j + 1
	}
	return( inner.products )
}
