get_g <- function(value, func_data){
    n<- length(value)
    j<- 0
    result <- rep(0,n)
    for(i in seq_len(n)){
        if(value[i]==func_data[j+1]){
            j = j+1
            result[i] = j - 1
        }else{
            result[i] = j - 1
        }
        if(j == length(func_data)&&i!=n){
            result[(i+1):n] = length(func_data)-1
            break
        }
    }
    result
}
get_h <- function(value, func_data){
    n<- length(value)
    j<- 0
    result <- rep(0,n)
    for(i in seq_len(n)){
        if(value[i]==func_data[j+1]){
            j = j+1
            result[i] = j
        }else{
            result[i] = j + 1 
        }
        if(j == length(func_data)&&i!=n){
            result[(i+1):n] = length(func_data)+1
            break
        }
    }
    result
}

#' @export
orderedProb1 <- function(l,h){
    if(length(l)==0) return(NA)
    if(any(l>=h)) return(0)
    n <- length(l)
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    
    prob <- KSgeneral:::compute_noncrossing_prob(h,l)
    return(prob)
}


#' @export
orderedProb <- function(l,h){
    if(length(l)==0) return(NA)
    if(any(l>=h)) return(0)
    n <- length(l)
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    
    total <- unique(sort(c(0,l,h,1)))
    #g(t_i)
    g_value <- get_g(total,h)
    #h(t_i)
    h_value <- get_h(total,l)
    
    n_t <- length(total)
    diff_t <- diff(total)
    m <- length(l)
    compute_prob_fft(m,g_value,h_value,n_t,diff_t)
}

orderedProb2 <- function(l,h){
    if(length(l)==0) return(NA)
    if(any(l>=h)) return(0)
    n <- length(l)
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    
    total <- unique(sort(c(0,l,h,1)))
    #g(t_i)
    g_value <- get_g(total,h)
    #h(t_i)
    h_value <- get_h(total,l)
    
    n_t <- length(total)
    diff_t <- diff(total)
    m <- length(l)
    compute_prob_fft2(m,g_value,h_value,n_t,diff_t)
}


increasingSeq <- function(from, to){
    if(to>=from){
        from:to
    }else{
        c()
    }
}
# orderedProb2(L,U)
# L <- c(0.2,0.5, 0.6,0.8);U <- c(0.4,0.7, 0.8,0.9)
# L <- c(0.2);U <- c(0.4)
# l = L;h=U

## h(L_i) = i
## g(U_i) = i - 1
## h(L_i + delta) = i + 1
## g(U_i + delta) = i - 1
## h(0) = 1
## h(1) = n + 1
## g(0) = - 1
## g(1) = n - 1
## S = T(g) U T(h) U {1}
orderedProbNative <- function(L,U){
    n <- length(L)
    T_g <- U
    T_h <- L
    S <- sort(unique(c(T_g,T_h,1)))
    ## Q[i] = Q(t, i - 1)
    Q <- rep(0, n+1)
    Q_new <- rep(0, n + 1)
    Q[1] <- 1
    t_i <- 0
    g_t_i <- -1
    # g_value <- get_g(S,T_g)
    #h(t_i)
    # h_value <- get_h(S,T_h)
    ## Fix t first
    ## then fix m
    for(i in seq_along(S)){
        t_i_plus_1 <- S[i]
        idx <- which(T_g<=t_i_plus_1)
        if(length(idx)){
            g_t_i_plus_1 <- max(which(T_g<=t_i_plus_1))-1
        }else{
            g_t_i_plus_1 <- -1
        }
        idx <- max(which(T_h<=t_i_plus_1))
        if(T_h[idx]==t_i_plus_1){
            h_t_i_plus_1 <- idx
        }else{
            h_t_i_plus_1 <- idx + 1
        }
        m_range <- increasingSeq(g_t_i_plus_1+1, h_t_i_plus_1-1)
        for(m in m_range){
            l_range <- increasingSeq(g_t_i + 1,m)
            Q_new[m + 1] <-  
                sum(Q[l_range + 1] * dpois(m-l_range, n*(t_i_plus_1-t_i)))
        }
        ## prepare for the next iteration
        # print(Q_new)
        Q <- Q_new
        Q_new <- rep(0, n + 1)
        t_i <- t_i_plus_1
        g_t_i <- g_t_i_plus_1
    }
    Q[n+1] / dpois(n,n)
}

