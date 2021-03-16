#include<stdint.h>
#include <vector>
#include <cmath>

std::vector<double> myFactorial;

double getPoisson(double k, double rate)
{
	if (rate == 0)
	{
		if (k == 0)
		{
			return (1);
		}
		else
		{
			return (0);
		}
	}
	double result = k * std::log(rate) - rate - myFactorial[k];
	return std::exp(result);
}

void computeFactorialUpTo(uint64_t n)
{
	uint64_t oldN = myFactorial.size();
	myFactorial.reserve(n + 1);
	for (uint64_t i = oldN; i < n + 1; ++i)
	{
		if (i == 0 || i == 1)
		{
			myFactorial[i] = 0;
		}
		else
		{
			myFactorial[i] = myFactorial[i - 1] + std::log(i);
		}
	}
}