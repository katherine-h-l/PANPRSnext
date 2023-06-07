#include <Rcpp.h>

// [[Rcpp::export]]
std::vector<int> FibCpp0(int n)
{
  // Error checking
  if(n <= 0)
    {
      throw std::range_error("n must be a positive integer");
    }

  // Allocate memory
  std::vector<int> out(n);
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
