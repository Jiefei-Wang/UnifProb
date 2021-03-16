#include <Rcpp.h>
#include <chrono>
#include "fftwconvolver.h"
#include <complex.h>
#include "fftw3.h"
#include "utils.h"
#include "Double_ptr.h"
#include <stdint.h>
using namespace Rcpp;

R_xlen_t findMaxMemSize(R_xlen_t n, double *g_value, double *h_value)
{
	R_xlen_t maxSize = 0;
	for (R_xlen_t i = 0; i < n - 1; ++i)
	{
		R_xlen_t curSize = h_value[i + 1] - g_value[i] - 1;
		maxSize = maxSize > curSize ? maxSize : curSize;
	}
	maxSize *= 2;
	double logSize = std::log2((double)maxSize);
	logSize = std::ceil(logSize);
	maxSize = std::pow(2, logSize);
	return maxSize;
}
R_xlen_t findCurMemSize(double g, double h)
{
	R_xlen_t curSize = h - g - 1;
	curSize *= 2;
	double logSize = std::log2((double)curSize);
	logSize = std::ceil(logSize);
	curSize = std::pow(2, logSize);
	return curSize;
}
void convolution(R_xlen_t N,
				 double *x_real, double *x_img,
				 double *y_real, double *y_img);

template <class T>
void convolution_noFFT(R_xlen_t N, T *out, T *buffer, T *x, T *y)
{
	for (R_xlen_t i = 0; i < N; ++i)
	{
		buffer[i] = 0;
		for (R_xlen_t j = 0; j <= i; ++j)
		{
			T curX = x[j];
			T curY = y[i - j];
			buffer[i] += curX * curY;
		}
	}
	memcpy(out, buffer, N * sizeof(T));
}

#define MINIMUM_CONVOLUTION_SIZE 80
// [[Rcpp::export]]
double compute_prob_fft(R_xlen_t m, SEXP R_g_value, SEXP R_h_value,
						R_xlen_t n_t, SEXP R_diff_t)
{
	computeFactorialUpTo(m);
	double *g_value = (double *)DATAPTR(R_g_value);
	double *h_value = (double *)DATAPTR(R_h_value);
	double *diff_t = (double *)DATAPTR(R_diff_t);

	R_xlen_t max_size = findMaxMemSize(n_t, g_value, h_value);
	//double* x_real = new double[max_size];
	double *Q = new double[m + 1 + max_size];
	double *y_real = new double[max_size];
	double *Q_img = new double[max_size];
	double *y_img = new double[max_size];

	SETVALUE(Q, m + max_size + 1, 0);
	Q[0] = 1;
	for (R_xlen_t i = 0; i < n_t - 1; ++i)
	{
		double gt_i = g_value[i];
		double ht_i_plus = h_value[i + 1];
		R_xlen_t maxRange = ht_i_plus - gt_i - 1;
		if (maxRange >= MINIMUM_CONVOLUTION_SIZE)
		{
			R_xlen_t curSize = findCurMemSize(gt_i, ht_i_plus);
			//Rprintf("Curlen:%d,MaxMem:%d\n", maxRange, curSize);
			SETVALUE(y_real + maxRange, curSize - maxRange, 0);
			SETVALUE(Q_img, curSize, 0);
			SETVALUE(y_img, curSize, 0);

			for (R_xlen_t j = 0; j < maxRange; j++)
			{
				//x_real[j] = Q[j + (R_xlen_t)gt_i + 1];
				y_real[j] = getPoisson(j, m * diff_t[i]);
			}
			convolution(curSize, Q + (R_xlen_t)gt_i + 1, Q_img, y_real, y_img);
			SETVALUE(Q + (R_xlen_t)gt_i + 1 + maxRange, curSize - maxRange, 0);
		}
		else
		{
			for (R_xlen_t j = 0; j < maxRange; j++)
			{
				y_real[j] = getPoisson(j, m * diff_t[i]);
			}
			convolution_noFFT(maxRange, Q + (R_xlen_t)gt_i + 1, Q_img, Q + (R_xlen_t)gt_i + 1, y_real);
		}
		//Rprintf("\n");
	}
	double result = Q[m] / getPoisson(m, m);
	delete[] Q;
	delete[] Q_img;
	delete[] y_real;
	delete[] y_img;
	return result;
}

// [[Rcpp::export]]
double compute_prob_fft2(R_xlen_t m, NumericVector &g_value, NumericVector &h_value,
						 R_xlen_t n_t, NumericVector &diff_t)
{
	computeFactorialUpTo(m);
	R_xlen_t max_size = findMaxMemSize(n_t, g_value.begin(), h_value.begin());

	FFTWConvolver fftconvolver(max_size);
	//double* x_real = new double[max_size];
	//NumericVector Q(m+1 + max_size);
	Double_ptr Q(m + 1 + max_size);
	//Buffers for doing the fft
	double *poisson_buffer = (double *)fftw_malloc(max_size * sizeof(double));
	double *Q_buffer = (double *)fftw_malloc(max_size * sizeof(double));

	//Initialize Q
	SETVALUE(Q.get_ptr(), Q.length, 0);
	SETVALUE(Q.get_another_ptr(), Q.length, 0);
	Q.get_ptr()[0] = 1;
	for (R_xlen_t i = 0; i < n_t - 1; ++i)
	{
		double &gt_i = g_value[i];
		double &ht_i_plus = h_value[i + 1];
		R_xlen_t maxRange = ht_i_plus - gt_i - 1;
		R_xlen_t curSize = findCurMemSize(gt_i, ht_i_plus);
		size_t start_Q_offset = (R_xlen_t)gt_i + 1;
		//memcpy(Q_buffer, Q.begin()+start_Q_offset, maxRange*sizeof(double));
		for (R_xlen_t j = 0; j < maxRange; j++)
		{
			//Q_buffer[j] = Q[j + start_Q_offset];
			poisson_buffer[j] = getPoisson(j, m * diff_t[i]);
		}
		//convolution(curSize, Q + (R_xlen_t)gt_i + 1, Q_img, y_real, y_img);
		fftconvolver.convolve_same_size(curSize, Q.get_ptr() + start_Q_offset, poisson_buffer, Q.get_another_ptr() + start_Q_offset);
		//convolution(curSize, Q + (R_xlen_t)gt_i + 1, Q_img, y_real, y_img);
		SETVALUE(Q.get_another_ptr() + start_Q_offset + maxRange, curSize - maxRange, 0);
		Q.switch_ptr();
		//Rprintf("\n");
	}
	double result = Q.get_another_ptr()[m] / getPoisson(m, m);
	fftw_free(poisson_buffer);
	fftw_free(Q_buffer);
	return result;
}

// [[Rcpp::export]]
NumericVector simpleConvolve(NumericVector &input1, NumericVector &input2)
{
	R_xlen_t n = Rf_length(input1);
	double *input_buffer1 = (double *)fftw_malloc(n * sizeof(double));
	double *input_buffer2 = (double *)fftw_malloc(n * sizeof(double));
	double *output_buffer = (double *)fftw_malloc(n * sizeof(double));

	memcpy(input_buffer1, DATAPTR(input1), n * sizeof(double));
	memcpy(input_buffer2, DATAPTR(input2), n * sizeof(double));
	FFTWConvolver convolver(n);
	convolver.convolve_same_size(n, input_buffer1, input_buffer2, output_buffer);

	NumericVector output(n);
	memcpy(DATAPTR(output), output_buffer, n * sizeof(double));

	fftw_free(input_buffer1);
	fftw_free(input_buffer2);
	fftw_free(output_buffer);
	return output;
}

#include "Convolver.h"

uint64_t fft_rounding = 1;
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
double compute_prob_fft3(R_xlen_t m, NumericVector &gt, NumericVector &ht,
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
			for (size_t k = 0; k < required_len; k++)
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