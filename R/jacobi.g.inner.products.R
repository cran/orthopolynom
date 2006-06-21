jacobi.g.inner.products <- function( n, p, q )
{
###
###   This function returns a vector with n+1 elements
###   containing the inner product of an order k Jacobi polynomial
###   Gk(p,q,x) with itself (i.e. norm squared) for orders k=0,1,...,n
###
###   Parameters
###   n = integer highest polynomial order
###   p = first parameter
###   q = second parameter
###
   if ( n < 0 )
      stop( "negative highest polynomial order" )
   if ( n != round( n ) )
      stop( "highest polynomial order is not integer" )
   if ( ( p - q ) <= -1 )
      stop( "p minus q less than or equal to -1" )
   if ( q <= 0 )
      stop( "q less than or equal to 0" )
   inner.products <- rep( 1, n + 1 )
   pq <- p - q
   j <- 1
   for ( k in 0:n ) {
      log.num <- lfactorial( k ) + lgamma( k + q ) + lgamma( k + p ) + lgamma( k + pq + 1 )
      log.den <- 2 * lgamma( 2 * k + p )
      inner.products[j] <- exp( log.num - log.den ) / ( 2 * k + p )
      j <- j + 1
   }
   return( inner.products )
}
