library(matlib)

funkcijasvd = function(A) {
  m = nrow(A)  #broj vrsta
  n = ncol(A)  #broj kolona
  
  #A transponovano puta A
  AtA = t(A) %*% A
  #dekompozicija matrice A trans puta A, preko eigen 
  eig = eigen(AtA)
  
  #eigenvalues su sopstvene vrednosti i smestaju se u lambda
  #sigma singularne vrednosti 
  lambda = eig$values
  lambda[lambda<0] = 0 #problem pri pokretanju Nans produced, pa treba ovaj uslov
  sigma = sqrt(lambda)
  
  #gram schmidt za ortonormirani skup 
  V = GramSchmidt(eig$vectors, normalize = TRUE)
  
  #rang matrice broj sing vrednosti vecih od nula, prag je 1e-8 da se i te jako male tretiraju kao nula 
  tol = 1e-8
  r = sum(sigma > tol)
  
 #prazna matrica kasnije ce biti popunenja
  U = matrix(0, m, m)
  
  #Računa leve singularne vektore (kolone matrice U). Samo za prvih r kolona jer su ostale nula, a ne moze se deliti s nulom
  for (i in 1:r) {
    U[, i] = (A %*% V[, i]) / sigma[i]
  }
  
  # ako je r manje od m onda matrica nema dovoljno kolona
  if (r < m) {
    kandidati = cbind(U[, 1:r, drop = FALSE], diag(m)) #dodau se jedinice
    Q = GramSchmidt(kandidati, normalize = TRUE) #pa ih gram schmidt ortogonalizuje
    U = Q[, 1:m] #matrica U je mxm sa svim kolonama ortonomiranim
  }
  
  
  #matrica D matrica singularnih vrednosti 
  D = matrix(0, m, n)
  for (i in 1:min(m, n)) {
    D[i, i] = sigma[i] #na dijagonali sing. vrednostsi 
  }
  
  sigmavec = sigma[1:min(m, n)] #izbacuje se eventualni visak 
  
  list(
    d = sigmavec ,
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

rezultat = funkcijasvd(A)


zapsmall(rezultat$d)
zapsmall(svd(A)$d)
razlika = rezultat$d - svd(A)$d
print(razlika) #razlika izmedju ove fje i svd su brojevi na e-14 e-15 e-8 e-33 e-34 vrlo bliski nuli tj mala je razlika
Arekonstr = rezultat$u %*% rezultat$D %*% t(rezultat$v)
print(round(Arekonstr))
print(A)


