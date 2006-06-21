hermite.he.polynomials <- function( n, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k scaled Hermite polynomials, He-k(x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### normalized = a boolean value.  if true, the polynomials are normalized
###
    recurrences <- hermite.he.recurrences( n, normalized )
    polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
    return( NULL )
}
