
# two sided
set.seed(1)
n <- 4000
x <-runif(n)
indexL <- seq_len(n)
indexU <- seq_len(n)
BJ <- UnifProb:::BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- UnifProb:::BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

## KS general
system.time(
    for(i in 1:10)
    orderedProb1(bounds$l,bounds$h)
)
## our method native
system.time(
    for(i in 1:10)
        UnifProb:::orderedProb(bounds$l,bounds$h)
)
## our method with FFTW
system.time(
    for(i in 1:10)
        UnifProb:::orderedProb2(bounds$l,bounds$h)
)

# n = 4000
# KS general
# user  system elapsed 
# 1.14    0.00    1.21 
# our method native
# user  system elapsed 
# 15.04    0.02   15.34 
# our method with FFTW
# user  system elapsed 
# 2.28    0.00    2.34 


# one sided
set.seed(1)
n <- 4000
x <-runif(n)
indexL <- seq_len(n)
indexU <- NULL
BJ <- BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

system.time(
    for(i in 1:2)
        orderedProb(bounds$l,bounds$h)
)


system.time(
    for(i in 1:2)
        orderedProb1(bounds$l,bounds$h)
)


# n = 4000
# reference
# user  system elapsed 
# 1.98    0.00    2.00 
# version 1
# user  system elapsed 
# 46.44    0.02   46.78 
# version 2
# user  system elapsed 
# 46.44    0.02   23.42 
