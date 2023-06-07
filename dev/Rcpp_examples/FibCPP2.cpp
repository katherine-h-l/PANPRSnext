// This example is identical to FibCpp1 with
// the exception that we use a namespace.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector FibCpp2(int n)
{
  // Error checking
  if(n <= 0)
    {
      throw std::range_error("n must be a positive integer");
    }

  // Allocate memory
  IntegerVector out(n);
  out[0]=1;

  // Compute additional terms
  if(n > 0)
    {
      out[1]=1;
      int i;
      for(i=2; i<n; i++)
	{
	  out[i] = out[i-1] + out[i-2];
	}
    }

  return out;
}
