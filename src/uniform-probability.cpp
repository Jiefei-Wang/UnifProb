#include <Rcpp.h>
#include <chrono>
#include <complex.h>
#include <stdint.h>
#include "fftw3.h"
#include "utils.h"
#include "Convolver.h"
using namespace Rcpp;

uint64_t fft_rounding = 256;
uint64_t fft_min_size = 80;

uint64_t find_fft_size(uint64_t n)
{
	bool remainder = (2 * n - 1) % fft_rounding;
	uint64_t size = ((2 * n - 1) / fft_rounding + remainder) * fft_rounding;
	return size;
}
uint64_t find_max_range(NumericVector &gt, NumericVector &ht)
{
	uint64_t n = gt.length();
	uint64_t max_range = 0;
	for (uint64_t i = 0; i < n - 1; i++)
	{
		max_range = max_range > (ht[i + 1] - gt[i]) ? max_range : (ht[i + 1] - gt[i]);
	}
	return max_range + 1;
}

// [[Rcpp::export]]
double compute_prob_fft(R_xlen_t m, NumericVector &gt, NumericVector &ht,
						 NumericVector &diff_t, bool debug = false)
{
	computeFactorialUpTo(m);
	uint64_t n_t = diff_t.length();
	uint64_t max_range = find_max_range(gt, ht);
	uint64_t max_fft_range = find_fft_size(max_range);

	/*
	steps:
	Q values -> buffer_in -> r2c -> buffer_out1
	poisson values -> buffer_in -> r2c -> buffer_out2
	buffer_out1 * buffer_out2 / length -> buffer_out1
	buffer_out1 -> c2r -> buffer_in
	*/
	double *buffer_in = fftw_alloc_real(max_fft_range);
	fftw_complex *buffer_out1 = fftw_alloc_complex(max_fft_range / 2 + 1);
	fftw_complex *buffer_out2 = fftw_alloc_complex(max_fft_range / 2 + 1);
	//Initialize Q
	buffer_in[0] = 1;
	SETVALUE(buffer_in + 1, max_fft_range - 1, 0);

	//This is for the non-fft algorithm
	NumericVector tmp_poisson(max_range);
	NumericVector tmp_out(max_range);

	Convolver convolver(buffer_in, buffer_out1, buffer_out2);
	for (uint64_t i = 0; i < n_t; ++i)
	{
		uint64_t cur_max_range = ht[i + 1] - gt[i] + 1;
		uint64_t required_len = ht[i + 1] - gt[i + 1] + 1;
		if (cur_max_range > fft_min_size)
		{
			uint64_t fft_n = find_fft_size(cur_max_range);

			convolver.set_size(fft_n);
			//Do the FFT on Q
			convolver.fft_out1();
			//Fill the poisson values
			for (uint64_t j = 0; j < cur_max_range; j++)
			{
				buffer_in[j] = getPoisson(j, m * diff_t[i]);
			}
			//Do the FFT on the poisson values
			convolver.fft_out2();
			//Convolution
			convolver.convolve();
			//Get the Q values for the next iteration
			uint64_t offset = gt[i + 1] - gt[i];
			if (offset != 0)
			{
				memmove(buffer_in, buffer_in + offset, required_len * sizeof(double));
			}
			//Set the rest to 0
			SETVALUE(buffer_in + required_len, fft_n - required_len, 0);
		}
		else
		{
			for (uint64_t j = 0; j < cur_max_range; j++)
			{
				tmp_poisson[j] = getPoisson(j, m * diff_t[i]);
			}
			for (uint64_t n = 0; n < cur_max_range; n++)
			{
				tmp_out[n] = 0;
				for (uint64_t k = 0; k <= n; k++)
				{
					tmp_out[n] += buffer_in[k]*tmp_poisson[n-k];
				}
			}
			uint64_t offset = gt[i + 1] - gt[i];
			memcpy(buffer_in, tmp_out.begin() + offset, required_len * sizeof(double));
			SETVALUE(buffer_in + required_len, cur_max_range - required_len, 0);
		}
		if (debug)
		{
			Rprintf("Iteration %d: ", i);
			for (size_t k = 0; k < 6; k++)
			{
				Rprintf("%f,", buffer_in[k]);
			}
			Rprintf("\n");
		}
	}
	double result = buffer_in[0] / getPoisson(m, m);
	fftw_free(buffer_in);
	fftw_free(buffer_out1);
	fftw_free(buffer_out2);
	return result;
}

#include "fftw_plan_manager.h"
// [[Rcpp::export]]
void set_plan_flag(Rcpp::String flag)
{
	if (flag == "FFTW_ESTIMATE")
	{
		fftw_plan_manager::set_flag(FFTW_ESTIMATE);
		return;
	}
	if (flag == "FFTW_MEASURE")
	{
		fftw_plan_manager::set_flag(FFTW_MEASURE);
		return;
	}
	if (flag == "FFTW_PATIENT")
	{
		fftw_plan_manager::set_flag(FFTW_PATIENT);
		return;
	}
	if (flag == "FFTW_EXHAUSTIVE")
	{
		fftw_plan_manager::set_flag(FFTW_EXHAUSTIVE);
		return;
	}
	if (flag == "FFTW_WISDOM_ONLY")
	{
		fftw_plan_manager::set_flag(FFTW_WISDOM_ONLY);
		return;
	}
	Rf_error("The flag does not exist\n");
}
// [[Rcpp::export]]
void set_fft_rounding(uint64_t x)
{
	fft_rounding = x;
}
// [[Rcpp::export]]
void set_fft_min_size(uint64_t x)
{
	fft_min_size = x;
}