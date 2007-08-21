jacobi.g.polynomials <- function( n, p, q, normalized=FALSE )
{
###
###   This function returns a list with n+1 elements
###   containing the order k Jacobi polynomials Gk(p,q,x)
###   for orders k=0,1,...,n
###
###   Parameters
###   n = integer highest polynomial order
###   p = first polynomial parameter
###   q = second polynomial parameter
###   normalized = boolean value.  if true, the polynomials are normalized
###
    require( polynom )
    almost.slegendre <- ( abs( p - q ) ) < 1e-6 & ( abs ( q - 1 ) < 1e-6 )
    if ( almost.slegendre )
        return( slegendre.polynomials( n, normalized ) )
    if ( (p != 2 ) && ( q != 1 ) ) {
       recurrences <- jacobi.g.recurrences( n, p, q, normalized )
       polynomials <- orthogonal.polynomials( recurrences )
    }
    else {
        np1 <- n + 1
        polynomials <- as.list( rep( NULL, np1 ) )
        k <- 0
        kp1 <- 1
        while ( k <= n ) {
            d.k <- exp( lgamma( k + 1 ) - lgamma( 2 * k + 2 ) )
            coef.k <- rep( 0, kp1 )
            sign.m <- 1
            m <- 0
            while ( m <= k ) {
                c.mmk <- sign.m * choose( k, m ) * exp( lgamma( 2 + 2 * k - m ) - lgamma( 1 + k - m ) )
                coef.k[kp1-m] <- c.mmk
                sign.m <- (-1) * sign.m
                m <- m + 1
            }
            polynomials[[kp1]] <- polynomial( coef.k )
            k <- k + 1
            kp1 <- kp1 + 1
        }
        if ( normalized ) {
            norms <- sqrt( jacobi.g.inner.products( n , p, q ) )
            for ( k in 1:np1 ) {
                polynomials[[k]] <- polynomials[[k]] / norms[k]
            }
        }
    }
    return( polynomials )
}
