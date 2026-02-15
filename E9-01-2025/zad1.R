library(matlib)

svd_manual <- function(A) {
  
  # dimenzije matrice
  m <- nrow(A)
  n <- ncol(A)
  
  #racunamo A^T A
  AtA <- t(A) %*% A
  
  #eigen dekompozicija
  eig <- eigen(AtA)
  
  #singularne vrednosti
  lambda <- eig$values
  lambda[lambda < 0] <- 0       # numericka stabilnost
  sigma <- sqrt(lambda)
  
  #desni singularni vektori (V)
  V <- GramSchmidt(eig$vectors, normalize = TRUE)
  
  #rang matrice
  tol <- 1e-8
  r <- sum(sigma > tol)
  
  #leva singularna matrica U
  U <- matrix(0, m, m)
  
  for (i in 1:r) {
    U[, i] <- (A %*% V[, i]) / sigma[i]
  }
  
  #nema dovoljno kolona → dopuna baze
  if (r < m) {
    kandidati <- cbind(U[, 1:r, drop = FALSE], diag(m))
    U <- GramSchmidt(kandidati, normalize = TRUE)[, 1:m]
  }
  
  #dijagonalna matrica singularnih vrednosti
  D <- matrix(0, m, n)
  for (i in 1:min(m, n)) {
    D[i, i] <- sigma[i]
  }
  
  list(
    d = sigma[1:min(m, n)],
    u = U,
    v = V,
    D = D
  )
}

A = matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12
), ncol = 6, byrow = TRUE)

result = svd_manual(A)


zapsmall(result$d)
zapsmall(svd(A)$d)
difference = result$d - svd(A)$d
print(difference) #razlika izmedju fje i svd su brojevi vrlo bliski nuli
Arekonstr = result$u %*% result$D %*% t(result$v)
print(round(Arekonstr))
print(A)
