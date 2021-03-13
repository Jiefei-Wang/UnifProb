#include <Rcpp.h>
#include <vector>
#include <chrono>
#include "fftwconvolver.h"
#include "fftw3.h"
using namespace Rcpp;

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

void computeFactorialUpTo(R_xlen_t n)
{
	R_xlen_t oldN = myFactorial.size();
	myFactorial.reserve(n + 1);
	for (R_xlen_t i = oldN; i < n + 1; ++i)
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

#define SETVALUE(x, n, value)        \
	for (R_xlen_t I = 0; I < (R_xlen_t)n; ++I) \
	{                                \
		(x)[I] = value;              \
	}

template <class T>
SEXP allocWithInit(int type, R_xlen_t len, T value)
{
	SEXP R_res = Rf_allocVector(type, len);
	T *res = (double *)DATAPTR(R_res);
	for (R_xlen_t i = 0; i < len; ++i)
	{
		res[i] = value;
	}
	return R_res;
}

// [[Rcpp::export]]
double compute_prob(R_xlen_t m, SEXP R_g_value, SEXP R_h_value,
					R_xlen_t n_t, SEXP R_diff_t)
{
	computeFactorialUpTo(m);
	double *Q = new double[m + 1];
	SETVALUE(Q, m + 1, 0);
	Q[0] = 1;
	double *newQ = new double[m + 1];
	double *g_value = (double *)DATAPTR(R_g_value);
	double *h_value = (double *)DATAPTR(R_h_value);
	double *diff_t = (double *)DATAPTR(R_diff_t);
	for (R_xlen_t i = 0; i < n_t - 1; ++i)
	{
		double gt_i = g_value[i];
		double gt_i_plus = g_value[i + 1];
		double ht_i_plus = h_value[i + 1];
		R_xlen_t maxRange = ht_i_plus - gt_i_plus - 1;
		SETVALUE(newQ, maxRange, 0);
		for (R_xlen_t j = 0; j < maxRange; ++j)
		{
			R_xlen_t curM = j + gt_i_plus + 1;
			for (R_xlen_t l = gt_i + 1; l <= curM; ++l)
			{
				double curQ = Q[l];
				double curPi = getPoisson(curM - l, m * diff_t[i]);
				newQ[j] += curQ * curPi;
			}
		}
		for (R_xlen_t j = gt_i_plus + 1; j <= ht_i_plus - 1; ++j)
		{
			Q[j] = newQ[j - (R_xlen_t)gt_i_plus - 1];
		}
		/*for(int k=0;k<m+1;k++){
			Rprintf("%f,", Q[k]);
		}
		Rprintf("\n");*/
	}
	double result = Q[m] / getPoisson(m, m);
	delete[] Q;
	delete[] newQ;
	return result;
}

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

#include "Double_ptr.h"
#include "fftwconvolver.h"
// [[Rcpp::export]]
double compute_prob_fft2(R_xlen_t m, NumericVector &g_value, NumericVector &h_value,
						 R_xlen_t n_t, NumericVector &diff_t)
{
	computeFactorialUpTo(m);
	R_xlen_t max_size = findMaxMemSize(n_t, g_value.begin(), h_value.begin());

	FFTWConvolver fftconvolver(max_size);
	//double* x_real = new double[max_size];
	//NumericVector Q(m+1 + max_size);
	Double_ptr Q(m+1 + max_size);
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
		fftconvolver.convolve_same_size(curSize, Q.get_ptr()+start_Q_offset, poisson_buffer, Q.get_another_ptr()+ start_Q_offset);
		//convolution(curSize, Q + (R_xlen_t)gt_i + 1, Q_img, y_real, y_img);
		SETVALUE(Q.get_another_ptr()+ start_Q_offset + maxRange, curSize - maxRange, 0);
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