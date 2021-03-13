#include <fftw3.h>
#include "Double_ptr.h"
Double_ptr::Double_ptr(size_t n) : length(n), is_switched(false)
{
    ptr1 = (double *)fftw_malloc(length * sizeof(double));
    ptr2 = (double *)fftw_malloc(length * sizeof(double));
}
Double_ptr::~Double_ptr()
{
    fftw_free(ptr1);
    fftw_free(ptr2);
}

void Double_ptr::switch_ptr(){
    is_switched = !is_switched;
}

double* Double_ptr::get_ptr(){
    return is_switched?ptr2:ptr1;
}
double* Double_ptr::get_another_ptr(){
    return is_switched?ptr1:ptr2;
}
