l <- c(0.1,0.1,0.7)
h <- c(0.3,0.4,0.9)
# set_fft_min_size(0)
orderedProbCFFT(l,h)


orderedProbNative(l,h)


set_fft_rounding(512)
a <- UnifProb:::orderedProbFFT(bounds$l,bounds$h,debug=F)
set_fft_rounding(128)
b <- UnifProb:::orderedProbFFT(bounds$l,bounds$h,debug=F)