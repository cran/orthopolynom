slegendre.polynomials <- function( n, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k shifted Legendre polynomials Pk(x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### normalized = boolean value.  if true, the polynomials are normalized
###
    recurrences <- slegendre.recurrences( n, normalized )
    polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
