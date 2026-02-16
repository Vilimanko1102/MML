library(matlib)

A = matrix(c(2, 1, 
             1, 20), nrow = 2, byrow = TRUE)
b = c(5, 3)

x_curr = c(-3, -1)      
gamma = 0.08         
n_iterations = 100    

for (i in 1:n_iterations) {
  grad_row = t(x_curr) %*% A - b
  x_curr = x_curr - gamma * t(grad_row)
}

print("Rešenje nakon 100 iteracija:")
print(x_curr)

print("Teorijsko rešenje (x* = A^-1 * b):")
print(solve(A) %*% b)

