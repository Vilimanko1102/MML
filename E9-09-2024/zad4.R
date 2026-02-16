library(matlib)

n = 4
x = rep(0.2, n)           
gamma = 0.01              
n_iter = 100         

A <- matrix(0, n, n)
for (i in 1:n) {
  for (j in i:n) {
    A[i, j] = 1
  }
}

f <- function(x) {
  Ax = A %*% x
  term1 = t(x) %*% (t(A) %*% A) %*% x
  term2 = 0.1 * sum(Ax^3)
  term3 = 0.01 * sum(Ax^4)
  return(term1[1] + term2 + term3)  
}

grad_f <- function(x) {
  Ax = A %*% x
  grad1 = 2 * t(A) %*% A %*% x
  grad2 = 0.3 * t(A) %*% (Ax^2)
  grad3 = 0.04 * t(A) %*% (Ax^3)
  return(grad1 + grad2 + grad3)
}

for (k in 1:n_iter) {
  grad = grad_f(x)
  x = x - gamma * as.vector(grad)
}

sprintf("Spall-ov polinom - rezultati nakon %d iteracija", n_iter)
print("x* = ") 
print(round(x, 6))
sprintf("f(x*) = %f", f(x))
