gegenbauer.recurrences <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a data frame with n+1 and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k generalized Laguerre polynomials, Lk(alpha,x),
### and for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
### normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( alpha <= -0.5 )
        stop( "alpha less than or equal to -0.5" )
    if ( alpha == 0 )
        return( chebyshev.t.recurrences( n, normalized ) )
    if ( alpha == 0.5 )
        return( legendre.recurrences( n, normalized ) )
    if ( alpha == 1 )
        return( chebyshev.u.recurrences( n, normalized ) )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    if ( normalized ) {
        norms <- sqrt( gegenbauer.inner.products( np1, alpha ) )
        while ( j <= n ) {
            r[k,"c"] <- ( j + 1 ) * norms[k+1]
            r[k,"d"] <- 0
            r[k,"e"] <- 2 * ( j + alpha ) * norms[k]
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                if ( k == 1 )
                    r[k,"f"] <- 0
                else
                    r[k,"f"] <- ( j + 2 * alpha - 1 ) * norms[k-1]
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            r[k,"c"] <- j + 1
            r[k,"d"] <- 0
            r[k,"e"] <- 2 * ( j + alpha )
            r[k,"f"] <- j + 2 * alpha - 1
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}
