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
