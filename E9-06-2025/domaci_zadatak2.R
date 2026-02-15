#Kolokvijum 1 — resenje koristeci my_svd()

my_svd <- function(A) {
  A <- as.matrix(A)
  
  ATA <- t(A) %*% A
  eg <- eigen(ATA, symmetric = TRUE)
  
  ord <- order(eg$values, decreasing = TRUE)
  lambda <- eg$values[ord]
  V <- eg$vectors[, ord, drop = FALSE]
  
  lambda[lambda < 0] <- 0
  d <- sqrt(lambda)
  
  r <- sum(d > 0)
  
  Vr <- V[, 1:r, drop = FALSE]
  dr <- d[1:r]
  
  U <- A %*% Vr
  U <- sweep(U, 2, dr, "/")
  
  list(d = dr, u = U, v = Vr)
}

#---------------------------------------------
#ZADATAK 1

A <- matrix(c( 2,  3,  0,  1,
               -1,  1,  2,  0,
               0,  2,  1,  1,
               1, -1, -4,  2,
               3,  0, -2, -1),
            nrow = 5,
            byrow = TRUE)
y <- matrix(c(-8, 3, -5, -17, 1), ncol = 1)

cat("ZADATAK 1\n")
cat("(a):\n")
#(a)
sA <- my_svd(A)
dA <- sA$d
rA <- length(dA)

cat("rank(A) = ", rA, "\n")
cat("dim ImΦ = ", rA, "\n")
cat("dim KerΦ = ", ncol(A) - rA, "\n\n")

#Baza ImΦ: kolone 1,2,3
A0 <- A[, 1:3, drop = FALSE]
cat("Baza ImΦ:\n")
print(A0)

#Baza KerΦ:
ATA <- t(A) %*% A
eg <- eigen(ATA, symmetric = TRUE)

ord <- order(eg$values, decreasing = TRUE)
lambda <- eg$values[ord]
Vfull <- eg$vectors[, ord, drop = FALSE]

lambda[lambda < 0] <- 0
d_full <- sqrt(lambda)
r_full <- sum(d_full > 0)

Ker_basis <- Vfull[, (r_full+1):ncol(A), drop = FALSE]
cat("\nBaza KerΦ:\n")
print(Ker_basis)

cat("\n(b):")
#(b) Koordinate y u bazi ImΦ
sA0 <- my_svd(A0)
U0 <- sA0$u
V0 <- sA0$v
d0 <- sA0$d
r0 <- length(d0)

A0_plus <- V0 %*% diag(1/d0, nrow = r0, ncol = r0) %*% t(U0)
c_coords <- A0_plus %*% y

cat("\nKoordinate y u bazi ImΦ (A0-baza):\n")
print(c_coords)

cat("\nProvera A0*c - y:\n")
print(round(A0 %*% c_coords - y, 10))

#(c) Jedno resenje x1
cat("\n(c):")
x1 <- rbind(c_coords, 0)
cat("\nJedno resenje x1:\n")
print(x1)

cat("\nProvera A*x1 - y :\n")
print(round(A %*% x1 - y, 10))

#(c) Drugo resenje:
if (ncol(Ker_basis) > 0) {
  k <- Ker_basis[,1, drop = FALSE]
  x2 <- x1 + k
  cat("\nDrugo resenje x2 = x1 + k:\n")
  print(x2)
  
  cat("\nProvera A*x2 - y:\n")
  print(round(A %*% x2 - y, 10))
}

#-------------------------------------------
#ZADATAK 2

cat("\n\nZADATAK 2\n")

M<-matrix(c( 0, -1, 2, 0, 2,
            -3, 4, 1, 2, 1,
            1, -3, 1,- 1, 2,
            -1, -3, 5, 0, 7),ncol=4); M
x <- matrix(c(-9,-1,-1, 4, 1), ncol = 1)


#SVD(M)
sM <- my_svd(M)
UM <- sM$u
dM <- sM$d
rM <- length(dM)

cat("dim U = rank(M) = ", rM, "\n")

proj <- UM %*% (t(UM) %*% x)
cat("\n(a) Projekcija pi_U(x):\n")
print(proj)

dist <- sqrt(sum((x - proj)^2))
cat("\n(b) Udaljenost d(x,U):\n")
print(dist)

#--------------------------------
#ZADATAK 3

cat("\n\nZADATAK 3\n")

A3 <- matrix(c(4,-3,-2,
               2,-1,-2,
               3,-3,-1), nrow = 3, byrow = TRUE)

v1 <- matrix(c(1,0,1), ncol = 1)
v2 <- matrix(c(1,1,0), ncol = 1)
v3 <- matrix(c(1,1,1), ncol = 1)

cat("A v1 =\n"); print(A3 %*% v1)
cat("A v2 =\n"); print(A3 %*% v2)
cat("A v3 =\n"); print(A3 %*% v3)

lambda1 <- as.numeric((t(v1) %*% (A3 %*% v1)) / (t(v1) %*% v1))
lambda2 <- as.numeric((t(v2) %*% (A3 %*% v2)) / (t(v2) %*% v2))
lambda3 <- as.numeric((t(v3) %*% (A3 %*% v3)) / (t(v3) %*% v3))

cat("\nlambda1=",lambda1," lambda2=",lambda2," lambda3=",lambda3,"\n")

P <- cbind(v1,v2,v3)
D <- diag(c(lambda1, lambda2, lambda3))
Pinv <- solve(P)

cat("\nProvera A - P D P^{-1}:\n")
print(A3 - P %*% D %*% Pinv)
