jacobi.p.recurrences <- function( n, a, b, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coeffieicnts c, d, e and f of the recurrence relations
### for the order k Jacobi polynomial, Pk(a,b,x),
### and for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### a = first polynomial parameter
### b = second polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( a <= -1 )
        stop( "alpha less than or equal to -1" )
    if ( b <= -1 )
        stop( "beta less than or equal to -1" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    ab <- a + b
    aabb <- a * a - b * b
    if ( normalized ) {
        norms <- sqrt( jacobi.p.inner.products( n+1, a, b ) )
        while ( j <= n ) {
            c <- 2 * ( j + 1 ) * ( j + ab + 1 ) * ( 2 * j + ab )
            d <- ( 2 * j + ab + 1 ) * aabb
            e <- pochhammer( 2 * j + ab, 3 )
            f <- 2 * ( j + a ) * ( j + b ) * ( 2 * j + ab + 2 )
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
            c <- 2 * ( j + 1 ) * ( j + ab + 1 ) * ( 2 * j + ab )
            d <- ( 2 * j + ab + 1 ) * aabb
            e <- pochhammer( 2 * j + ab, 3 )
            f <- 2 * ( j + a ) * ( j + b ) * ( 2 * j + ab + 2 )
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
