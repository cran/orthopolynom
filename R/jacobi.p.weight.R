jacobi.p.weight <- function( x, alpha, beta )
{
###
###	This function returns the value of the weight function
###	for the Jacobi polynomial, Pk( alpha, beta, x )
###
###	Parameters
###	x = the function argument
###	alpha = the first polynomial parameter
###	beta = the second polynomial parameter
###
	n <- length( x )
	y <- rep( 0, n )
	for ( i in 1:n ) {
		if ( ( x[i] > -1 ) && ( x[i] < +1 ) )
			y[i] <- ( ( 1 - x[i] ) ^ alpha ) * ( ( 1 + x[i] ) ^ beta )
	}
	return ( y )
}