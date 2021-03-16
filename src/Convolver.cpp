#include <complex>
#include "Convolver.h"
#include "fftw_plan_manager.h"


Convolver::Convolver(double *buffer_in, fftw_complex *buffer_out1, fftw_complex *buffer_out2) : buffer_in(buffer_in), buffer_out1(buffer_out1), buffer_out2(buffer_out2) {}

void Convolver::set_size(uint64_t n)
{
    size = n;
}
void Convolver::fft_out1()
{
    fftw_plan_manager::exec_r2c(size, buffer_in, buffer_out1);
}
void Convolver::fft_out2()
{
    fftw_plan_manager::exec_r2c(size, buffer_in, buffer_out2);
}
void Convolver::convolve()
{
    std::complex<double> *buffer1 = reinterpret_cast<std::complex<double> *>(buffer_out1);
    std::complex<double> *buffer2 = reinterpret_cast<std::complex<double> *>(buffer_out2);
    double constant = 1 / (double)size;
    for (uint64_t i = 0; i < size / 2 + 1; i++)
    {
        buffer1[i] = buffer1[i] * buffer2[i] * constant;
    }
    fftw_plan_manager::exec_c2r(size, buffer_out1, buffer_in);
}