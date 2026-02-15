library(matlib)
#prekopirana iz zad 1.
funkcijasvd = function(A) {
  m = nrow(A)
  n = ncol(A)
  AtA = t(A) %*% A
  eig = eigen(AtA)
  lambda = eig$values
  lambda[lambda<0] = 0 
  sigma = sqrt(lambda)
  V = GramSchmidt(eig$vectors, normalize = TRUE)
  tol = 1e-8
  r = sum(sigma > tol)
  U = matrix(0, m, m)
  for (i in 1:r) {
    U[, i] = (A %*% V[, i]) / sigma[i]
  }
  if (r < m) {
    kandidati = cbind(U[, 1:r, drop = FALSE], diag(m)) 
    Q = GramSchmidt(kandidati, normalize = TRUE) 
    U = Q[, 1:m]
  }
  D = matrix(0, m, n)
  for (i in 1:min(m, n)) {
    D[i, i] = sigma[i] 
  }
  sigmavec = sigma[1:min(m, n)] 
  list(
    d = sigmavec,
    u = U,
    v = V,
    D = D
  )
}


# 1. zadatak 
#zbog sukoba istih prom koriscene su oznake sa brojem zadatke


z1_A = matrix(c(2,3,0,1, -1,1,2,0, 0,2,1,1, 1,-1,-4,2, 3,0,-2,-1), ncol=4, byrow=T)
z1_y = matrix(c(-8, 3, -5, -17, 1))

z1_svd = funkcijasvd(z1_A)
z1_r = sum(z1_svd$d > 1e-8)

#dim i baza im i ker
cat("dim Im Phi", z1_r, "\n")
cat("dim Ker Phi", ncol(z1_A) - z1_r, "\n")

z1_a0 = z1_A[, 1:3];
z1_a0 #bazaim

z1_ker = matrix(c(-1,1,-1,-1)) #baza ker
z1_ker

# izraziti y u bazi Im Phi
z1_svdA0 = funkcijasvd(z1_a0)
z1_r0 = sum(z1_svdA0$d > 1e-8)

z1_Dplus0 = matrix(0, ncol(z1_a0), nrow(z1_a0))
for(i in 1:z1_r0) {
  z1_Dplus0[i,i] = 1/z1_svdA0$d[i]
}

#pseudoinverz
z1_a10 = z1_svdA0$v %*% z1_Dplus0 %*% t(z1_svdA0$u)

#koeficijenti y u bazi Im Phi
z1_x10 = z1_a10 %*% z1_y; 
z1_x10 
# 5 * prva kolona baze -6 puta druga kolona baze + 7 puta treca kolona baze 

#provera da y pripada Im Phi
round(z1_a0 %*% z1_x10 - z1_y, digits = 8)  #treba sve nule da se dobiju 


z1_x1 = rbind(z1_x10, 0); 
z1_x1  # x1 = (5, -6, 7, 0)

round(z1_A %*% z1_x1 - z1_y, digits = 8) 

#jedinstveno? Ne jer dim Ker = 1
z1_x0 = matrix(c(-1, 1, -1, -1)); 
z1_x0 

z1_x2 = z1_x1 + z1_x0;
z1_x2  # x2 = (4, -5, 6, -1)

round(z1_A %*% z1_x2 - z1_y, digits = 8)  # Provera

#x3 nekolinearno sa x2 - x1 ne postoji jer dim Ker = 1



#2. zadatak


z2_U = matrix(c(0,1,-3,-1,-1,-3,4,-3,2,1,1,5,0,-1,2,0,2,2,1,7), ncol=4, byrow=TRUE);
z2_x = matrix(c(-9,-1,-1,4,1), ncol=1);
z2_x;

#projekcija x na U preko svd pseudoinverza
z2_svd = funkcijasvd(z2_U)
z2_rU = sum(z2_svd$d > 1e-8)

z2_Dplus = matrix(0, ncol(z2_U), nrow(z2_U))
for(i in 1:z2_rU) {
  z2_Dplus[i,i] = 1/z2_svd$d[i]
}

z2_Uplus = z2_svd$v %*% z2_Dplus %*% t(z2_svd$u)

z2_Px = z2_U %*% z2_Uplus %*% z2_x;
z2_Px;

#rastojanje
z2_xPx = matrix(z2_x - z2_Px, ncol=1);

sqrt(t(z2_xPx) %*% z2_xPx);  # 6.728287


#3. zadatak 


#
z3_A = matrix(c(4,-3,-2,2,-1,-2,3,-3,-1), ncol=3, byrow=TRUE);

z3_x1 = matrix(c(1,0,1), ncol=1);
z3_x2 = matrix(c(1,1,0), ncol=1);
z3_x3 = matrix(c(1,1,1), ncol=1);

z3_l1 = 2;
z3_l2 = 1;
z3_l3 = -1;

# Provera: A*v - lambda*v treba da je nula 
z3_A %*% z3_x1 - z3_l1 * z3_x1;  # treba sve nule, x1 jeste kar. vektor za lambda=2
z3_A %*% z3_x2 - z3_l2 * z3_x2;  # treba sve nule, x2 jeste kar. vektor za lambda=1
z3_A %*% z3_x3 - z3_l3 * z3_x3;  # treba sve nule, x3 jeste kar. vektor za lambda=-1

# rekosntrukcija A za proveru
z3_svd = funkcijasvd(z3_A)
round(z3_svd$u %*% z3_svd$D %*% t(z3_svd$v), 6)  # rekonstrukcija A

#dijagonalizacija A = P * D * P^(-1)
z3_P = cbind(z3_x1, z3_x2, z3_x3);
z3_D = diag(c(z3_l1, z3_l2, z3_l3));
z3_D;
z3_P;

# preko svd pseudoioznerva
z3_svdP = funkcijasvd(z3_P)

z3_DplusP = matrix(0, 3, 3)
for(i in 1:3) {
  z3_DplusP[i,i] = 1/z3_svdP$d[i]
}

z3_P1 = z3_svdP$v %*% z3_DplusP %*% t(z3_svdP$u);
z3_P1;

#provera
z3_provera = z3_P %*% z3_D %*% z3_P1;
z3_provera;

round((z3_provera - z3_A),8) # treba sve nule