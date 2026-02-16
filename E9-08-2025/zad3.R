 

#f(x) = 0.5 [x1, x2] · [2, 1; 1, 20] · [x1, x2] - [5, 3] · [x1, x2]
#                             A                     b

#matrica A i vektor b iz primera se direktno procitaju
A = matrix(c(2, 1, 1, 20), ncol = 2, byrow = TRUE)
b = c(5, 3)

#pocetna tacka
 x = c(-3, -1)

#korak 
gamma = 0.07

#broj iteracija
brojiter = 100

#tacke za ispis 
putanja = matrix(0, brojiter + 1, 2)
putanja[1, ] = x

#gradijentni postupak
for (i in 1:brojiter) {
  grad = A %*% x - b           #gradjient
  x = x - gamma * grad         #nova tacka
  putanja[i + 1, ] = x
}

# rezultati
cat("Pocetna tacka: ", putanja[1, ], "\n")
cat("Krajnja tacka posle 100 iteracija: ", putanja[brojiter + 1, ], "\n")



cat("\n Svih 100 iteracija:\n")
for (i in 1:100) {
  cat(round(putanja[i, 2], 6),"\n")
}

