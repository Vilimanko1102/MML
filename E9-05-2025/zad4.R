# Problem 5.1)

n = 4
A = matrix(c(1,0,0,0,
             1,1,0,0,
             1,1,1,0,
             1,1,1,1), ncol = n)

f = function(x) {
  as.numeric(t(x) %*% t(A) %*% A %*% x + 0.1 * sum(A %*% x^3) + 0.01 * sum(A %*% x^4))
}

grad_f = function(x) {
  2 * t(A) %*% A %*% x + 0.3 * t(A) %*% (A %*% x^2) + 0.04 * t(A) %*% (A %*% x^3)
}

x = matrix(c(0.2, 0.2, 0.2, 0.2), ncol = 1)
alpha = 0.1

hist = data.frame(i = integer(), x1 = numeric(), x2 = numeric(), x3 = numeric(), x4 = numeric(), fx = numeric())
for (i in 1:100) {
  grad = grad_f(x)
  x = x - alpha * grad
  if (i %% 10 == 0) {
    hist = rbind(hist, data.frame(i = i, x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], fx = f(x)))
  }
}

print(hist, digits = 13) # x -> [0,0,0,0]; f(x) -> 0
