#Domaci zadatak 4
# Problem 5.1 – Spall's Polynomial
# Primeniti Gradient Descent (100 iteracija)

#dimenzije
n <- 4

#gornja trougaona matrica A
A <- matrix(0, n, n)
A[upper.tri(A, diag = TRUE)] <- 1

#ciljna fja
f <- function(x) {
  y <- A %*% x
  sum(y^2) + 0.1 * sum(y^3) + 0.01 * sum(y^4)
}

#gradijent fje
grad_f <- function(x) {
  y <- A %*% x
  g_y <- 2*y + 0.3*y^2 + 0.04*y^3
  t(A) %*% g_y
}

#gradijentni spust
gradient_descent <- function(x0, gamma, n_iter) {
  x <- x0
  history <- matrix(NA, nrow = n_iter + 1, ncol = n)
  history[1, ] <- x
  
  for (k in 1:n_iter) {
    x <- x - gamma * grad_f(x)
    history[k + 1, ] <- x
  }
  
  history
}

#pocetna tacka x0
x0 <- 0.2 * rep(1, n)

#parametri learning rate i broj ponavljanja
gamma <- 0.01
n_iter <- 100

#uradi gradijentni spust
trajectory <- gradient_descent(x0, gamma, n_iter)

trajectory[101, ]
