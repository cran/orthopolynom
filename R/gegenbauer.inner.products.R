gegenbauer.inner.products <- function( n, alpha )
{
###
### This function returns a vector with n+1 elements
### containing the inner product of an order k Gegenbauer polynomial, Ck(alpha,x),
### with itself (i.e. norm squared) for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
###
    if ( n < 0 )
        stop( "highest polynomial order is less than zero" )
    if ( n != round(n) )
        stop( "highest polynomial order is not an integer" )
    if ( alpha <= -0.5 )
        stop( "alpha is less than or equal to -0.5" )
    if ( abs( alpha ) < 1e-6 )
        return( chebyshev.t.inner.products( n ) )
    if ( abs( alpha - 0.5 ) < 1e-6 )
        return( legendre.inner.products( n ) )
    if ( abs( alpha - 1.0 ) < 1e-6 )
        return( chebyshev.u.inner.products( n ) )
    inner.products <- rep( 1, n + 1 )
    j <- 1
    coef <- pi * ( 2 ^ ( 1 - 2 * alpha ) )
    for ( k in 0:n ) {
        if ( abs( alpha ) < 1e-6 ) {
            if ( k == 0 )
                inner.products[j] <- pi
            else
                inner.products[j] <- ( 2 * pi ) / ( k ^ 2 )
        }
        else {
            log.num <- lgamma( k + 2 * alpha )
            log.den <- lfactorial( k ) + 2 * lgamma( alpha )
            if ( k == -alpha ) {
                    inner.products[j] <- coef * exp( log.num - log.den )
            }
            else {
                inner.products[j] <- coef * exp( log.num - log.den ) / ( k + alpha )
            }
        }
        j <- j + 1
    }
    return( inner.products )
}
