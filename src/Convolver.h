#include <fftw3.h>
#include <stdint.h>
class Convolver
{
	double *buffer_in;
	fftw_complex *buffer_out1;
	fftw_complex *buffer_out2;
	uint64_t size;
public:
	Convolver(double *buffer_in, fftw_complex *buffer_out1, fftw_complex *buffer_out2);

	void set_size(uint64_t n);
	void fft_out1();
	void fft_out2();
	void convolve();
};