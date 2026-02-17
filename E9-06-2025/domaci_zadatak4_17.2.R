#1)parametri 

n <- 4

#gornja trougaona matrica popunjena jedinicama
A <- matrix(0, n, n)
for (i in 1:n) {
  for (j in i:n) {
    A[i, j] <- 1
  }
}

#----------------------------
#2)Funkcija cilja f(x)

f <- function(x) {
  x <- as.numeric(x)
  Ax <- A %*% x
  
  term1 <- t(x) %*% t(A) %*% A %*% x
  term2 <- 0.1 * sum(Ax^3)
  term3 <- 0.01 * sum(Ax^4)
  
  as.numeric(term1 + term2 + term3)
}

#----------------------------
#3)gradijent

grad_f <- function(x) {
  x <- as.numeric(x)
  Ax <- A %*% x
  
  term1 <- 2 * t(A) %*% A %*% x
  term2 <- 0.3 * t(A) %*% (Ax^2)
  term3 <- 0.04 * t(A) %*% (Ax^3)
  
  as.numeric(term1 + term2 + term3)
}

#----------------------------
#4)gradient descent


x <- rep(0.2, n)   #početna tačka
alpha <- 0.001     #learning rate

num_iter <- 100

#istorija za praćenje
history_x <- matrix(NA, nrow = num_iter + 1, ncol = n)
history_f <- numeric(num_iter + 1)

history_x[1, ] <- x
history_f[1] <- f(x)

for (k in 1:num_iter) {
  g <- grad_f(x)
  x <- x - alpha * g
  
  history_x[k + 1, ] <- x
  history_f[k + 1] <- f(x)
}

#----------------------------
#5) Ispis rezultata

cat("Rezultat posle", num_iter, "iteracija \n")
cat("x =", paste(round(x, 6), collapse = ", "), "\n")
cat("f(x) =", format(history_f[num_iter + 1], digits = 10), "\n\n")





