jacobi.g.weight <- function( x, p, q )
{
###
###   This function returns the value of the weight function
###   for the Jacobi polynomial, Gk( p, q, x )
###
###   Parameters
###   x = the function argument
###   p = the first polynomial parameter
###   q = the second polynomial parameter
###
    almost.slegendre <- ( abs( p - q ) ) < 1e-6 & ( abs ( q - 1 ) < 1e-6 )
    if ( almost.slegendre )
        return( slegendre.weight( x ) )
    n <- length( x )
    y <- rep( 0, n )
    for ( i in 1:n ) {
       if ( ( x[i] > 0 ) && ( x[i] < 1 ) ) {
          t1 <- ( p - q ) * log( 1 - x[i] )
          t2 <- ( q - 1 ) * log( x[i] )
          y[i] <- exp( t1 + t2 )
       }
    }
    return( y )
}
