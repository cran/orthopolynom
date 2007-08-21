jacobi.g.recurrences <- function( n, p, q, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k Jacobi polynomial, Gk(p,q,x),
### and for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### p = first polynomial parameter
### q = second polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    almost.slegendre <- ( abs( p - q ) ) < 1e-6 & ( abs ( q - 1 ) < 1e-6 )
    if ( almost.slegendre )
        return( slegendre.recurrences( n, normalized ) )
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( ( p - q ) <= -1 )
        stop( "p minus q less than or equal to -1" )
    if ( q <= 0 )
        stop( "q less than or equal to 0" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    if ( normalized ) {
        norms <- sqrt( jacobi.g.inner.products( n+1, p, q ) )
        while ( j <= n ) {
            c <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
            d <- -( 2 * j * ( j + p ) + q * ( p - 1 ) ) * pochhammer(2 * j + p - 2, 3 )
            e <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
            f <- j * ( j + q - 1 ) * ( j + p - 1 ) * ( j + p - q ) * ( 2 * j + p + 1 )
            r[k,"c"] <- c * norms[k+1]
            r[k,"d"] <- d * norms[k]
            r[k,"e"] <- e * norms[k]
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                if ( k == 1 )
                    r[k,"f"] <- 0
                else
                    r[k,"f"] <- f * norms[k-1]
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            c <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
            d <- -( 2 * j * ( j + p ) + q * ( p - 1 ) ) * pochhammer(2 * j + p - 2, 3 )
            e <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
            f <- j * ( j + q - 1 ) * ( j + p - 1 ) * ( j + p - q ) * ( 2 * j + p + 1 )
            r[k,"c"] <- c
            r[k,"d"] <- d
            r[k,"e"] <- e
            r[k,"f"] <- f
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}
