// This example is identical to FibCpp0 with
// the exception that std::vectors have been
// replaced by Rcpp::IntegerVector.

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector FibCpp1(int n)
{
  // Error checking
  if(n <= 0)
    {
      throw std::range_error("n must be a positive integer");
    }

  // Allocate memory
  Rcpp::IntegerVector out(n);
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
