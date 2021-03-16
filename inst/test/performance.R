
# two sided
set.seed(1)
n <- 8000
x <-runif(n)
indexL <- seq_len(n)
indexU <- seq_len(n)
BJ <- UnifProb:::BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- UnifProb:::BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

nsim <- 1
## ksgeneral
system.time(
    for(i in seq_len(nsim))
        UnifProb:::ksgeneral(bounds$l,bounds$h)
)

## our method with FFTW
# UnifProb:::set_plan_flag("FFTW_MEASURE")
system.time(
    for(i in seq_len(nsim))
        UnifProb:::orderedProbFFT(bounds$l,bounds$h)
)

# n = 4000
# our method native
# user  system elapsed 
# 15.04    0.02   15.34 
# our method with ksgeneral C++ code
# user  system elapsed 
# 2.28    0.00    2.34 
# KS general
# user  system elapsed 
# 1.14    0.00    1.21 


# one sided
set.seed(1)
n <- 4000
x <-runif(n)
indexL <- seq_len(n)
indexU <- NULL
BJ <- UnifProb:::BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- UnifProb:::BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

nsim <- 10

## ksgeneral
system.time(
    for(i in seq_len(nsim))
        UnifProb:::ksgeneral(bounds$l,bounds$h)
)

## our method with FFTW
# UnifProb:::set_plan_flag("FFTW_MEASURE")
system.time(
    for(i in seq_len(nsim))
        UnifProb:::orderedProbFFT(bounds$l,bounds$h)
)


# n = 4000
# reference
# user  system elapsed 
# 1.98    0.00    9.84
# version 1
# user  system elapsed 
# 46.44    0.02   46.78 
# version 2
# user  system elapsed 
# 46.44    0.02   23.42 
# version 3
# user  system elapsed 
# 46.44    0.02   6.72



# one sided
set.seed(1)
n <- 10000
x <-runif(n)
indexL <- seq_len(n)
indexU <- NULL
BJ <- UnifProb:::BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- UnifProb:::BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)
UnifProb:::ksgeneral(bounds$l,bounds$h)
set_fft_rounding(128)
UnifProb:::orderedProbFFT(bounds$l,bounds$h,debug=F)


nsim <- 1

## ksgeneral
system.time(
    for(i in seq_len(nsim))
        UnifProb:::ksgeneral(bounds$l,bounds$h)
)

## our method with FFTW
# UnifProb:::set_plan_flag("FFTW_MEASURE")
system.time(
    for(i in seq_len(nsim))
        UnifProb:::orderedProbFFT(bounds$l,bounds$h)
)
# Reference
# user  system elapsed 
# 8.29    0.00    8.30 
# version 1
# user  system elapsed 
# 5.64    0.00    5.64
