l <- c(0.1,0.1,0.7)
h <- c(0.3,0.4,0.9)
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

diff_t <- diff(total)
m <- length(l)
compute_prob_fft3(m,g_value,h_value,diff_t)


orderedProbNative(l,h)