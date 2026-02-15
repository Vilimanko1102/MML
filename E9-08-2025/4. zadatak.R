#zad 5.2 str 13.


b = 1

f = function(x) {
  s1 = sqrt(x[1]^2 + b * x[2]^2)
  s2 = sqrt(10^2 + b * 10^2)
  200 * (exp(0.1 * s1) + exp(-0.1 * s1) - 2) / (exp(0.1 * s2) + exp(-0.1 * s2) - 2)
}

grad_f = function(x) {
  h = 1e-6
  g1 = (f(c(x[1] + h, x[2])) - f(c(x[1] - h, x[2]))) / (2 * h)
  g2 = (f(c(x[1], x[2] + h)) - f(c(x[1], x[2] - h))) / (2 * h)
  return(c(g1, g2))
}

x = c(0.5, 0.5)
gamma = 0.1
brojiter = 100

putanja = matrix(0, brojiter + 1, 2)
putanja[1, ] = x

for (i in 1:brojiter) {
  grad = grad_f(x)
  x = x - gamma * grad
  putanja[i + 1, ] = x
}

print(paste("Start:", putanja[1,1], putanja[1,2]))
print(paste("Kraj:", round(putanja[brojiter+1,1],6), round(putanja[brojiter+1,2],6)))
print(paste("f(kraj):", round(f(putanja[brojiter+1,]),6)))
#poklapa se sa resenjm datim u knjizi 

