D <- 1
A <- 1
Is <- c(4)
Ks <- c(7)
f <- c(1)
arrX <- array(rep(0, prod(Is)), Is)
arrY <- array(rep(0, prod(Ks)), Ks)

for (i1 in 1:Is[1]) {
    arrX[i1] <- i1
}
for (k1 in 1:Ks[1]) {
    arrY[k1] <- k1
}

ps <- list()
for (a in 1:A) {
    ds <- which(f == a)
    d <- ds[1]
    length_of_p_a <- dim(arrX)[d]
    ps[[a]] <- rep(0.01, length_of_p_a); ps[[a]][c(1, 3)] <- 1
    ps[[a]] <- ps[[a]] / sum(ps[[a]])
}
qs <- list()
for (a in 1:A) {
    ds <- which(f == a)
    d <- ds[1]
    length_of_q_a <- dim(arrY)[d]
    qs[[a]] <- rep(1, length_of_q_a); qs[[a]][c(2, 3)] <- 0
    qs[[a]] <- qs[[a]] / sum(qs[[a]])
}

X <- as.tensor(arrX)
Y <- as.tensor(arrY)

out <- OTT(X = X, Y = Y, f = f,
    ps=ps, qs=qs, num.sample=2,
    loss = function (x, y) {abs(x - y)},
    num.iter=2, epsilon=1e-10)

expect_equal(length(out), 1)
expect_equal(length(out$Ts), 2)
expect_equal(dim(out$Ts[[1]]), c(4,6))
expect_equal(dim(out$Ts[[2]]), c(5,7))
