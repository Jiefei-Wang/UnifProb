
# two sided
set.seed(1)
n <- 4000
x <-runif(n)
indexL <- seq_len(n)
indexU <- seq_len(n)
BJ <- BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

system.time(
    for(i in 1:10)
    orderedProb(bounds$l,bounds$h)
)


system.time(
    for(i in 1:10)
    orderedProb1(bounds$l,bounds$h)
)



# n = 4000
# reference
# user  system elapsed 
# 0.97    0.00    0.97
# version 1
# user  system elapsed 
# 14.14    0.04   14.35 




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
