#Gradijentni postupak(100 iteracija)

A <- matrix(c(2, 1,
              1, 20), nrow = 2, byrow = TRUE)
b <- c(5, 3)

f <- function(x) {
  #f(x) = 1/2 x^T A x - b^T x
  0.5 * as.numeric(t(x) %*% A %*% x) - sum(b * x)
}

grad_f <- function(x) {
  #gradijent funkcije f(x) = A x - b
  as.numeric(A %*% x - b)
}

#Podešavanja
gamma <- 0.085          
n_iter <- 100
x <- c(-3, -1)          

#Čuvanje putanje iteracija
X <- matrix(NA_real_, nrow = n_iter + 1, ncol = 2)
F <- numeric(n_iter + 1)

X[1, ] <- x
F[1] <- f(x)

for (i in 1:n_iter) {
  x <- x - gamma * grad_f(x)
  X[i + 1, ] <- x
  F[i + 1] <- f(x)
}

#Rezultati
cat("x_100 =", X[n_iter + 1, ], "\n")
cat("f(x_100) =", F[n_iter + 1], "\n")

#analitičko rešenje (pošto je funkcija kvadratna)
x_star <- solve(A, b)
cat("x* (analitičko) =", x_star, "\n")
cat("f(x*) =", f(x_star), "\n")

#prikaz prvih i poslednjih nekoliko iteracija
head(cbind(iter = 0:n_iter, X, f = F), 6)
tail(cbind(iter = 0:n_iter, X, f = F), 6)
