get_g_t<- function(h, x){
    h <- c(h,1.1)
    result <- rep(0, length(x))
    j <- 1
    for(i in seq_along(x)){
        while (x[i]>=h[j]) {
            j=j+1
        }
        result[i] <- j - 1.5
    }
    return(result)
}
get_h_t<- function(l, x){
    l <- c(l, 1)
    result <- rep(0, length(x))
    j <- 1
    for(i in seq_along(x)){
        while (x[i]>l[j]) {
            j=j+1
        }
        result[i] <- j - 0.5
    }
    return(result)
}

## S = T(g) U T(h) U {1,0}
orderedProbNative <- function(l,h){
    n <- length(l)
    S <- sort(unique(c(0,1,l,h)))
    gt <- ceiling (get_g_t(h, S))
    ht <- floor(get_h_t(l, S))
    Q <- rep(0, n + 1)
    Q[1] <- 1
    ## Fix t first
    ## then fix m
    for(i in seq_len(length(S)-1)){
        Q_new <- rep(0, n + 1)
        t_i <- S[i]
        t_i_plus_1 <- S[i+1]
        g_i <- gt[i]
        g_i_plus_1 <- gt[i+1]
        h_i <- ht[i]
        h_i_plus_1 <- ht[i+1]
        m_range <- increasingSeq(g_i_plus_1, h_i_plus_1)
        for(m in m_range){
            l_range <- increasingSeq(g_i,m)
            Q_new[m + 1] <-  
                sum(Q[l_range + 1] * dpois(m-l_range, n*(t_i_plus_1-t_i)))
        }
        ## prepare for the next iteration
        print(Q_new)
        Q <- Q_new
    }
    Q[n+1] / dpois(n,n)
}


myconvolution <- function(x,y, strong= TRUE){
    nx <- length(x)
    ny <- length(y)
    stopifnot(nx==ny)
    if(strong){
        true_n <- (nx+1)/2
        if(nx!=1){
        stopifnot(nx%%2==1)
        stopifnot(all(x[(true_n+1):nx]==0))
        stopifnot(all(y[(true_n+1):nx]==0))
        }
    }
    x_fft <- fft(x)
    y_fft <- fft(y)
    z <- x_fft*y_fft/nx
    Re(fft(z, TRUE))
}

orderedProbNativeFFT <- function(l,h){
    n <- length(l)
    S <- sort(unique(c(0,1,l,h)))
    gt <- ceiling (get_g_t(h, S))
    ht <- floor(get_h_t(l, S))
    max_range <- max(ht[2:length(ht)]-gt[1:(length(gt)-1)])+1
    Q <- rep(0, 3*max_range-2)
    offset <- max_range
    Q[offset + 0] <- 1
    ## Fix t first
    ## then fix m
    for(i in seq_len(length(S)-1)){
        Q_new <- rep(0, 3*max_range-2)
        t_i <- S[i]
        t_i_plus_1 <- S[i+1]
        g_i <- gt[i]
        g_i_plus_1 <- gt[i+1]
        h_i <- ht[i]
        h_i_plus_1 <- ht[i+1]
        m_range <- g_i:h_i_plus_1
        m_num <- length(m_range)
        print(m_num)
        pois <- c(
            dpois(m_range - g_i, n*(t_i_plus_1-t_i)),
            rep(0, m_num - 1)
            )
        
        fft_n <- m_num*2-1
        Q_new_off <- offset + g_i - g_i_plus_1
        Q_new[Q_new_off+seq_len(fft_n)-1] <- 
            myconvolution(
                Q[offset + seq_len(fft_n)-1],
                pois
                )
        Q_new[increasingSeq(Q_new_off+m_num, Q_new_off+fft_n-1)] <- 0
        ## prepare for the next iteration
        print(Q_new)
        Q <- Q_new
    }
    Q[offset] / dpois(n,n)
}

orderedProbCFFT <- function(l,h){
    m <- length(l)
    S <- sort(unique(c(0,1,l,h)))
    gt <- ceiling (get_g_t(h, S))
    ht <- floor(get_h_t(l, S))
    dt <- diff(S)
    compute_prob_fft3(m,gt,ht,dt,TRUE)
}


