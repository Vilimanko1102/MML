#Domaci 2
library(matlib)

svd_manual <- function(A) {
  
  m <- nrow(A)
  n <- ncol(A)
  
  AtA <- t(A) %*% A
  eig <- eigen(AtA)
  
  lambda <- eig$values
  lambda[lambda < 0] <- 0
  sigma <- sqrt(lambda)
  
  V <- GramSchmidt(eig$vectors, normalize = TRUE)
  
  tol <- 1e-8
  r <- sum(sigma > tol)
  
  U <- matrix(0, m, m)
  for (i in 1:r) {
    U[, i] <- (A %*% V[, i]) / sigma[i]
  }
  
  if (r < m) {
    kandidati <- cbind(U[, 1:r, drop = FALSE], diag(m))
    U <- GramSchmidt(kandidati, normalize = TRUE)[, 1:m]
  }
  
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

#Zadatak 1.

A1 <- matrix(c(
  2,  3,  0,  1,
  -1,  1,  2,  0,
  0,  2,  1,  1,
  1, -1, -4,  2,
  3,  0, -2, -1
), ncol = 4, byrow = TRUE)

y1 <- matrix(c(-8, 3, -5, -17, 1))

svd_A1 <- svd_manual(A1)
r1 <- sum(svd_A1$d > 1e-8)

cat("dim Im Φ =", r1, "\n")
cat("dim Ker Φ =", ncol(A1) - r1, "\n")

A1_0 <- A1[, 1:3]
ker_A1 <- matrix(c(-1, 1, -1, -1))

svd_A1_0 <- svd_manual(A1_0)
r0 <- sum(svd_A1_0$d > 1e-8)

Dplus1 <- matrix(0, ncol(A1_0), nrow(A1_0))
for (i in 1:r0) {
  Dplus1[i, i] <- 1 / svd_A1_0$d[i]
}

A1_0_plus <- svd_A1_0$v %*% Dplus1 %*% t(svd_A1_0$u)
x_part <- A1_0_plus %*% y1

x1 <- rbind(x_part, 0)
x2 <- x1 + ker_A1

#Zadatak 2.

U2 <- matrix(c(
  0,  1, -3, -1,
  -1, -3,  4, -3,
  2,  1,  1,  5,
  0, -1,  2,  0,
  2,  2,  1,  7
), ncol = 4, byrow = TRUE)

x2_vec <- matrix(c(-9, -1, -1, 4, 1))

svd_U2 <- svd_manual(U2)
rU <- sum(svd_U2$d > 1e-8)

Dplus2 <- matrix(0, ncol(U2), nrow(U2))
for (i in 1:rU) {
  Dplus2[i, i] <- 1 / svd_U2$d[i]
}

U2_plus <- svd_U2$v %*% Dplus2 %*% t(svd_U2$u)
Px <- U2 %*% U2_plus %*% x2_vec
dist <- sqrt(t(x2_vec - Px) %*% (x2_vec - Px))

#Zadatak 3.

A3 <- matrix(c(
  4, -3, -2,
  2, -1, -2,
  3, -3, -1
), ncol = 3, byrow = TRUE)

v1 <- c(1, 0, 1); λ1 <- 2
v2 <- c(1, 1, 0); λ2 <- 1
v3 <- c(1, 1, 1); λ3 <- -1

P <- cbind(v1, v2, v3)
D <- diag(c(λ1, λ2, λ3))

svd_P <- svd_manual(P)
DplusP <- diag(1 / svd_P$d)
P_inv <- svd_P$v %*% DplusP %*% t(svd_P$u)

A_rec <- P %*% D %*% P_inv

round(A1 %*% x1 - y1, 8)
round(A1 %*% x2 - y1, 8)
round(A_rec - A3, 8)
dist
