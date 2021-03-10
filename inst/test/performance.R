BJPlusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    vapply(index, function(i)
        pbeta(sx[i], i, n - i + 1),numeric(1))
}
BJMinusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    1 - BJPlusLevel(n, x, sx, index)
}
BJStatFunc <- function(x, indexL=seq_along(x), indexU=seq_along(x)) {
    n <- length(x)
    sx <- sort(x)
    BJPlus <- min(BJPlusLevel(n, x, sx,indexL),1)
    BJMinus <- min(BJMinusLevel(n, x, sx,indexU),1)
    min(BJPlus, BJMinus)
}
BJLocalCritical<-function(statValue,n, indexL, indexU){
    l=vapply(seq_len(n),function(x)qbeta(statValue,x,n-x+1),numeric(1))
    h=vapply(seq_len(n),function(x)qbeta(1 - statValue,x,n-x+1),numeric(1))
    if(length(indexL)!=0){
        l[-indexL] <- 0
    }else{
        l=rep(0,length(l))
    }
    if(length(indexU)!=0){
        h[-indexU] <- 1
    }else{
        h=rep(1,length(h))
    }
    
    list(l =l,h= h)
}

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
