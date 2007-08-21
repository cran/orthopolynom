gegenbauer.polynomials <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Gegenbauer polynomial Ck(alpha,x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
### normalized = boolean value.  if true, the polynomials are normalized
###
    require( polynom )
    if ( abs( alpha ) < 1e-6 ) {
        t.polynomials <- chebyshev.t.polynomials( n, normalized )
        np1 <- n + 1
        polynomials <- as.list( rep( NULL, np1 ) )
        polynomials[[1]] <- 2 * t.polynomials[[1]]
        k <- 1
        kp1 <- 2
        while ( k <= n ) {
            t.polynomial <- t.polynomials[[kp1]]
            polynomials[[kp1]] <- ( 2 / k ) * t.polynomial
            k <- k + 1
            kp1 <- kp1 + 1
        }
    }
    else if ( abs( alpha - 0.5 ) < 1e-6 )
        polynomials <- legendre.polynomials( n, normalized )
    else if ( abs( alpha - 1.0 ) < 1e-6 )
        polynomials <- chebyshev.u.polynomials( n, normalized )
    else {
        recurrences <- gegenbauer.recurrences( n, alpha, normalized )
        if ( normalized ) {
            alpha.2 <- 2 * alpha
            h.0 <- pi * ( 2^(1-alpha.2) ) * exp ( lgamma(alpha.2) - lgamma(alpha) -lgamma(alpha+1) )
            p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
            polynomials <- orthonormal.polynomials( recurrences, p.0 )
        }
        else
            polynomials <- orthogonal.polynomials( recurrences )
    }
    return( polynomials )
}
