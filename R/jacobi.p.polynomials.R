jacobi.p.polynomials <- function( n, a, b, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Jacobi polynomials Pk(a,b,x)
### for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### a = first polynomial parameter
### b = second polynomial parameter
### normalized = boolean value.  if true, the polynomials are normalized
###
    require( polynom )
    almost.legendre <- ( abs( a ) < 1e-6 ) & ( abs( b ) < 1e-6 )
    if ( almost.legendre )
        return( legendre.polynomials( n, normalized ) )
    recurrences <- jacobi.p.recurrences( n, a, b, normalized )
    if ( normalized ) {
        ap1 <- a + 1
        bp1 <- b + 1
        abp1 <- a + b + 1
        h.0 <- ( ( 2^abp1 ) / abp1 ) * exp( lgamma(ap1 ) + lgamma(bp1) - lgamma(abp1) )
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
