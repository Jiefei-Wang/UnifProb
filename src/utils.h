#ifndef HEADER_UTILS
#define HEADER_UTILS
#include <Rcpp.h>
#include <stdint.h>


#define SETVALUE(x, n, value)                  \
	for (uint64_t macro_i = 0; macro_i < ((uint64_t)n); ++macro_i) \
	{                                          \
		(x)[macro_i] = value;                        \
	}



class PROTECT_GUARD
{
private:
  int protect_num = 0;

public:
  PROTECT_GUARD() {}
  ~PROTECT_GUARD()
  {
    if (protect_num != 0)
      UNPROTECT(protect_num);
  }
  SEXP protect(SEXP x)
  {
    protect_num++;
    return PROTECT(x);
  }
};



double getPoisson(double k, double rate);
void computeFactorialUpTo(uint64_t n);


#endif