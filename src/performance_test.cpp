#include <Rcpp.h>
#include<fftw3.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector performance_test1(NumericVector &input)
{
	R_xlen_t n = Rf_length(input);
	double *real_buffer_in = fftw_alloc_real(n);
	fftw_complex *complex_buffer = fftw_alloc_complex(n / 2 + 1);
	double *real_buffer_out_hidden = fftw_alloc_real(n);
	double *real_buffer_out = real_buffer_out_hidden;

	fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(n, real_buffer_in, complex_buffer, FFTW_ESTIMATE);
	fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(n, complex_buffer, real_buffer_out, FFTW_ESTIMATE);

	memcpy(real_buffer_in, DATAPTR(input), n * sizeof(double));
	fftw_execute(plan_r2c);
	fftw_execute(plan_c2r);
	NumericVector output(n);
	memcpy(DATAPTR(output), real_buffer_out, n * sizeof(double));

	fftw_free(real_buffer_in);
	fftw_free(complex_buffer);
	fftw_free(real_buffer_out_hidden);
	return output;
}

// [[Rcpp::export]]
NumericVector performance_test2(NumericVector &input)
{
	R_xlen_t n = Rf_length(input);
	double *real_buffer_in = fftw_alloc_real(n);
	fftw_complex *complex_buffer = fftw_alloc_complex(n / 2 + 1);
	double *real_buffer_out_hidden = fftw_alloc_real(n + 1);
	double *real_buffer_out = real_buffer_out_hidden + 1;

	fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(n, real_buffer_in, complex_buffer, FFTW_ESTIMATE);
	fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(n, complex_buffer, real_buffer_out, FFTW_ESTIMATE);

	memcpy(real_buffer_in, DATAPTR(input), n * sizeof(double));
	fftw_execute(plan_r2c);
	fftw_execute(plan_c2r);
	NumericVector output(n);
	memcpy(DATAPTR(output), real_buffer_out, n * sizeof(double));

	fftw_free(real_buffer_in);
	fftw_free(complex_buffer);
	fftw_free(real_buffer_out_hidden);
	return output;
}

// [[Rcpp::export]]
NumericVector performance_test3(NumericVector &input)
{
	R_xlen_t n = Rf_length(input);
	double *real_buffer_in_hidden = fftw_alloc_real(n + 1);
	fftw_complex *complex_buffer_hidden = fftw_alloc_complex(n / 2 + 1);
	double *real_buffer_out_hidden = fftw_alloc_real(n + 1);

	double *real_buffer_in = (double *)((char *)real_buffer_in_hidden + 1);
	fftw_complex *complex_buffer = (fftw_complex *)((char *)complex_buffer_hidden + 1);
	double *real_buffer_out = (double *)((char *)real_buffer_out_hidden + 1);

	fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(n, real_buffer_in, complex_buffer, FFTW_ESTIMATE);
	fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(n, complex_buffer, real_buffer_out, FFTW_ESTIMATE);

	memcpy(real_buffer_in, DATAPTR(input), n * sizeof(double));
	fftw_execute(plan_r2c);
	fftw_execute(plan_c2r);
	NumericVector output(n);
	memcpy(DATAPTR(output), real_buffer_out, n * sizeof(double));

	fftw_free(real_buffer_in_hidden);
	fftw_free(complex_buffer_hidden);
	fftw_free(real_buffer_out_hidden);
	return output;
}

// [[Rcpp::export]]
void performance_test4(R_xlen_t n)
{
	double *real_buffer_in = fftw_alloc_real(n);
	double *real_buffer_out = fftw_alloc_real(n);
	fftw_complex *complex_buffer = fftw_alloc_complex(n / 2 + 1);

	for (R_xlen_t i = 1; i <= n; i++)
	{
		fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(n, real_buffer_in, complex_buffer, FFTW_ESTIMATE);
		fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(n, complex_buffer, real_buffer_out, FFTW_ESTIMATE);
		fftw_destroy_plan(plan_r2c);
		fftw_destroy_plan(plan_c2r);
	}
	fftw_free(real_buffer_in);
	fftw_free(real_buffer_out);
	fftw_free(complex_buffer);
}