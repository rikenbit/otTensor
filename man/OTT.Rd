\name{OTT}
\docType{methods}
\alias{OTT}

\title{
Optimal Tensor Transport
}
\description{
Transport plans to align two tensors X and Y are estimated.
}
\usage{
OTT(X, Y, f, ps=NULL, qs=NULL,
    loss=.absolute_error, num.sample=1000,
    num.iter=200, epsilon=1e-10, verbose=FALSE)
}

\arguments{
  \item{X}{
The first tensor data ('rTensor' object). The order must be the same as that of Y.
}
  \item{Y}{
The second tensor data ('rTensor' object). The order must be the same as that of X.
}
  \item{f}{
Affectation function to assign transport plan to each pair of mode of X and Y.
}
  \item{ps}{
Row-wise weight vectors for transport plans (Default: NULL, which means uniform distribution).
}
  \item{qs}{
Column-wise weight vectors for transport plans (Default: NULL, which means uniform distribution).
}
  \item{loss}{
Loss function (Default: .absolute_error).
}
  \item{num.sample}{
Number of samples to calculate the gradient (Default: 1000).
}
  \item{num.iter}{
Number of iterations (Default: 200).
}
  \item{epsilon}{
Regularization parameter (Default: 1e-10).
}
  \item{verbose}{
Verbose option (Default: FALSE).
}
}

\value{
Ts : A list contains transport plans.
}
\author{Koki Tsuyuzaki}

\examples{
library("rTensor")
D <- 3
A <- 2
Is <- c(4, 4, 5)
Ks <- c(6, 6, 7)
f <- c(1, 1, 2)
arrX <- array(rep(0, prod(Is)), Is)
arrY <- array(rep(0, prod(Ks)), Ks)

for (i1 in 1:Is[1]) {
    for (i2 in 1:Is[2]) {
        for (i3 in 1:Is[3]) {
            arrX[i1, i2, i3] <- i1 + i2 + i3
        }
    }
}
for (k1 in 1:Ks[1]) {
    for (k2 in 1:Ks[2]) {
        for (k3 in 1:Ks[3]) {
            arrY[k1, k2, k3] <- k1 + k2 + k3
        }
    }
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

# Test Dataset
X <- as.tensor(arrX)
Y <- as.tensor(arrY)

# This is just for an example.
# In real data analysis,
# please specify larger num.sample and num.iter such as 1000 and 200, respectively.
OTT(X = X, Y = Y, f = f,
    ps=ps, qs=qs, num.sample=10,
    loss = function (x, y) {abs(x - y)},
    num.iter=2, epsilon=1e-10)
}

\keyword{methods}