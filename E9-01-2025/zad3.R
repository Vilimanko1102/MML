#Example 7.1.
#matrica A i vektor b
A <- matrix(c(2, 1,
              1, 20), nrow = 2, byrow = TRUE)

b <- c(5, 3)

#fja gradijenta
grad_f <- function(x) {
  A %*% x - b
}

#fja gradijentnog spusta
#x0  pocetna tacka, gamma brzina ucenja
gradient_descent <- function(x0, gamma, n_iter) {
  x <- x0
  history <- matrix(NA, nrow = n_iter + 1, ncol = 2)
  history[1, ] <- x
  
  for (k in 1:n_iter) {
    #racuna GD u trenutnoj tacki, pomeramo se suprotno od gradijenta
    x <- x - gamma * grad_f(x)
    history[k + 1, ] <- x #cuvamo novu poziciju
  }
  
  history
}

x0 <- c(-3, -1) #date vrednosti x0 iz knjige
gamma <- 0.085 #learning rate iz knjige
n_iter <- 100 #100 ponavljanja

trajectory <- gradient_descent(x0, gamma, n_iter)
trajectory[101,] #x100

x1 <- x0 - gamma * grad_f(x0) #proveravamo prvi korak
x1
solve(A, b) #GD kovergira do ove tacke

plot(trajectory[,1], trajectory[,2],
     type = "l", lwd = 2,
     xlab = "x1", ylab = "x2",
     main = "Gradient Descent Trajectory (Example 7.1)")

points(solve(A, b)[1], solve(A, b)[2],
       col = "red", pch = 19)
