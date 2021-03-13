test <- function(x,y, n_zeros){
    x_raw <- c(x,rep(0,n_zeros))
    y_raw <- c(y, rep(0,n_zeros))
    
    x_fft <- fft(x_raw)
    y_fft <- fft(y_raw)
    
    s <- x_fft*y_fft/length(x_raw)
    
    s1 <- Re(fft(s, TRUE))
    round(s1)
}




x <- c(1,2,3)
y <- c(4,5,6)

q_raw <- c(x,rep(0,length(x)))
p_raw <- c(y, rep(0,length(y)))

q_fft <- fft(q_raw)
p_fft <- fft(p_raw)

s <- q_fft*p_fft/length(q_raw)

s1 <- Re(fft(s, TRUE))

Re(s1)
s1[3]
x <- sum(Q[l_range + 1] * dpois(m-l_range, n*(t_i_plus_1-t_i)))

s1[3]/x


devtools::load_all()

x <- c(1,2,3)
y <- c(4,5,6)

simpleConvolve(x,y)




set.seed(1)
n <- 4000
x <-runif(n)
indexL <- seq_len(n)
indexU <- NULL
BJ <- BJStatFunc(x, indexL = indexL, indexU= indexU)
bounds <- BJLocalCritical(BJ, n, indexL = indexL, indexU= indexU)

system.time(
        a<-orderedProb(bounds$l,bounds$h)
)
system.time(
    b<-orderedProb2(bounds$l,bounds$h)
)
b<-orderedProb2(bounds$l,bounds$h)




compute_prob_fft2()