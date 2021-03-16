l <- c(0.1,0.1,0.7)
h <- c(0.3,0.4,0.9)

orderedProb(l,h)

orderedProb2(l,h)

orderedProbNative(l,h)




orderedProb1(l,h)

x <- c(1:3,0,0)
y <- c(4:6,0,0)
myconvolution(x,y,strong = FALSE)

x <- c(1:3,0,0,0)
y <- c(4:6,0,0,0)
myconvolution(x,y,strong = FALSE)

n <- 100000
x <- runif(10000000)
system.time({
    for(i in 1:10)
    res <- performance_test1(x)
})


system.time({
    for(i in 1:10)
    res2 <- performance_test2(x)
})

system.time({
    for(i in 1:10)
        res3 <- performance_test3(x)
})
res

system.time({
       performance_test4(10000)
})
