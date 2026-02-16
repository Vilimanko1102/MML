library(matlib)

# ===== ZADATAK 1 =====

A = matrix(c(2, -1, 0, 1, 3, 
             3,  1, 2, -1, 0, 
             0,  2, 1, -4, -2, 
             1,  0, 1, 2, -1), ncol = 4)
A

y = c(-8, 3, -5, -17, 1)

res = svd(A)
U = res$u
V = res$v
D = res$d

tol = 1e-10
rank_A = sum(D > tol)

# a)
basis_Im = U[, 1:rank_A]
basis_Im
dim_Im = rank_A
dim_Im
basis_Ker = V[, (rank_A + 1):ncol(A), drop = FALSE]
basis_Ker
dim_Ker = ncol(A) - rank_A
dim_Ker

# b)
y_coords_in_Im = t(basis_Im) %*% y
y_coords_in_Im

# c)
D_inv = diag(1/D[1:rank_A])
D_inv
x1 = V[, 1:rank_A] %*% D_inv %*% t(U[, 1:rank_A]) %*% y
x1

# ===== ZADATAK 2 =====

A2 = matrix(c( 0, -1,  2,  0,  2,
               1, -3,  1, -1,  2,
               -3,  4,  1,  2,  1,
               -1, -3,  5,  0,  7), 
            nrow = 5, ncol = 4)

x = c(-9, -1, -1, 4, 1)

res = svd(A2)
U_matrix = res$u
D = res$d

tol = 1e-10
rank_A = sum(D > tol)

basis_U = U_matrix[, 1:rank_A]

# a)
coords = t(basis_U) %*% x
pi_U_x = basis_U %*% coords
pi_U_x
# b)
residual_vector = x - pi_U_x
distance = sqrt(sum(residual_vector^2))
distance

# ===== ZADATAK 3 =====

A3 <- matrix(c(4, -3, -2, 
              2, -1, -2, 
              3, -3, -1), nrow = 3, byrow = TRUE)

P <- matrix(c(1, 0, 1,   # v1
              1, 1, 0,   # v2
              1, 1, 1),  # v3
            nrow = 3)

# a)
AP = A3 %*% P

lambdas = AP[1, ] / P[1, ]

is_eigen = all.equal(AP, P %*% diag(lambdas))
sprintf("Are v1, v2, v3 eigenvectors? %s", is_eigen)
sprintf("Eigenvalue : %s", unlist(lambdas))

# b)
D = diag(lambdas)
P1 = solve(P)
print("Diagonalization A = PDP^-1 is possible")
print("Matrix P:"); print(P)
print("Matrix D:"); print(D)
print("Matrix P1:"); print(P1)
