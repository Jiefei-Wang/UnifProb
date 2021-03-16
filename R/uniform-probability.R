orderedProbFFT <- function(l,h, debug = FALSE){
    m <- length(l)
    S <- sort(unique(c(0,1,l,h)))
    gt <- ceiling (get_g_t(h, S))
    ht <- floor(get_h_t(l, S))
    dt <- diff(S)
    compute_prob_fft(m,gt,ht,dt,debug)
}


ksgeneral <- function(l,h){
    if(length(l)==0) return(NA)
    if(any(l>=h)) return(0)
    n <- length(l)
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    
    prob <- KSgeneral:::compute_noncrossing_prob(h,l)
    prob
}

