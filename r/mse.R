library(Rssa)
source("toeplitz_mssa.R")

N <- 71
sigma <- 5
Ls <- c(12, 24, 36, 48, 60)
signal1 <- 30 * cos(2*pi * (1:N) / 12)
signal2 <- 30 * cos(2*pi * (1:N) / 12 + pi / 4)

signal <- cbind(signal1, signal2)
R <- 10000

ssa.errors <- function(LS) {
  f <- signal1 + rnorm(N, sd = sigma)
  err.rec <- numeric(length(Ls));
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "1d-ssa")
    rec <- reconstruct(s, groups = list(1:2))[[1]]
    err.rec[l] <- mean((rec - signal1)^2)
  }
  err.rec
}

toeplssa.errors <- function(Ls) {
  f <- signal1 + rnorm(N, sd = sigma)
  err.rec <- numeric(length(Ls));
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "toeplitz-ssa")
    rec <- reconstruct(s, groups = list(1:2))[[1]]
    err.rec[l] <- mean((rec - signal1)^2)
  }
  err.rec
}

mssa.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "mssa")
    rec <- reconstruct(s, groups = list(1:2))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
}

toeplSum.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- toeplitz.mssa(f, L = L, D = 2, method = "sum")
    rec <- toeplitz.reconstruct(s, groups = list(trend = c(1,2)))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
} 

toeplBlock.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- toeplitz.mssa(f, L = L, D = 2, method = "block")
    rec <- toeplitz.reconstruct(s, groups = list(trend = c(1,2)))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
} 

mres.ssa <- replicate(R, ssa.errors(Ls))
mres.toeplssa <- replicate(R, toeplssa.errors(LS))
err.rec <- rowMeans(mres)
err.rec





mse<-c()
for(i in 1:1000){
  f<-t + rnorm(N, sd=5)
    #cbind(h1 + rnorm(N,sd=5),h2+ rnorm(N,sd=5))
  s<-ssa(f,L = 60, kind = "toeplitz-ssa")
  r<-reconstruct(s,groups = list(trend=1:2))
  mse<-c(mse,mean((r$trend-t)^2))
}
mean(mse)









N <- 71
sigma <- 5
Ls <- c(12, 24, 36, 48, 60)
signal1 <- 6*(1:N)
signal2 <- 4*(1:N)
signal <- cbind(signal1, signal2)
R <- 10000
mssa.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- toeplitz.mssa(f, L = L, D = 2, method = "block")
    #s <- ssa(f, L = L, kind = "mssa")
    rec <- toeplitz.reconstruct(s, groups = list(trend = c(1,2)))$trend
    #rec <- reconstruct(s, groups = list(1:2))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
}
mres <- replicate(R, mssa.errors(Ls))
err.rec <- rowMeans(mres)
err.rec
