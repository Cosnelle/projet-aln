import numpy as np
import scipy.linalg as nla

# Affichage plus agréable
np.set_printoptions(linewidth=240)

A=np.array([[-7/10,-109/25,1/50],[7/10,-31/25,209/50],[23/10,-49/25,161/50 ],[-23/10,-91/25,49/50 ]])

m,n = A.shape

def reflecteur (p, v) :
  n = v.shape[0]
  F = np.eye (n) - 2 * np.outer (v,v)
  Q = np.eye (p, dtype=np.float64)
  Q [p-n:p,p-n:p] = F
  return Q

# Algorithme resolution Matrice triangulaire inférieure  substitution avant(forward)
def descente(A, b) :
    m = A.shape[0]
    for i in range (0, m) :
        for j in range (0, i) :
            b[i] = b[i] - A[i,j]*b[j]
        b[i] = b[i]/A[i,i]


# Algorithme resolution Matrice triangulaire suprieure  substitution avant(forward)
def remonte(A, b):
    m= A.shape[0]
    for i in range (m-1, -1, -1):
        for j in range (i+1, m):
            b[i] = b[i] - A[i,j]*b[j]
        b[i] = b[i]/A[i,i]



#E1
def bidiagonale(A):
   m,n = A.shape
   B = np.copy(A)
   VL = []
   VR = []
   for i in range (0,n):
        x = B[i:m,i]
        v = np.copy (x)
        v[0] = v[0]+np.sign(v[0])*nla.norm(x,2)
        v = (1/nla.norm(v,2))*v
        VL.append(v)
        Q = reflecteur (m, v)
        B = np.dot (Q,B)
        if(len(VR) < n-2):   
            x = B[i,i+1:n]
            v = np.copy(x)
            v[0] = v[0] + np.sign(v[0])*nla.norm(x,2)
            v = (1/nla.norm(v,2))*v
            VR.append(v)
            Q = reflecteur(n, v)
            B = np.dot (B, Q)
   B = B[0:n,0:n]
   return B, VL, VR

B, VL, VR = bidiagonale(A)
    

def calcul_svd(A):
m,n = A.shape
B, VL, VR = bidiagonale(A)
H = np.zeros([2*n, 2*n], dtype=np.float64)   #E2
H[0:n,n:2*n] = np.transpose(B)
H[n:2*n,0:n] = B
P = np.zeros ([2*n,2*n], dtype=np.float64)
for i in range (0,n) :
       P[i,2*i] = 1
       P[n+i,2*i+1] = 1
T = np.dot(np.transpose(P), np.dot(H,P))
d = np.zeros (2*n, dtype=np.float64)   #E3
e = np.array ([ T[i+1,i] for i in range (0,2*n-1) ], dtype=np.float64)
eigvals, eigvecs = nla.eigh_tridiagonal(d, e)
 
    #Algorithme QR pour trouver une approximation des valeurs propres
C = T
for k in range (0,100) :
      Q, R = nla.qr (C)
      C = np.dot (R, Q)

    #sur la diagonale de C on a une approximation des valeurs propres de T
def fi_signe(a,b,mu):
    n=len(a)
    f=[1]
    nbvp=0
    for i in range(1,(n+1)):
        alpha=a[i-1]-mu
        beta=b[i-1]*b[i-1]
        if b[i-1]==0 :
            f[i]=(alpha*np.sign(f[i-1]))
        elif b[i-2]==0 and b[i-1]!=0 :
            f[i]=(alpha*f[i-1]-beta*np.sign(f[i-2]))
        else :
            f[i]=(alpha*f[i-1]-beta*f[i-2])
        if f[i]!=0:
            f.append(f[i])
        else:
            f[i].append(f[i-1])
        if f[i]*f[i-1]<0:
            nbvp=nbvp+1
    return f,nbvp

m= T.shape[0]
a = []
for i in range (0,m):
  a.append(T[i,i])

b= [0]
for i in range(0,m-1):
      b.append(T[i+1, i])


m = len(b)
max = 0
for i in range(0,m):
    if(abs(b[i]) > max):
        max = b[i]

maxbi = abs(max)

fi_signe(a,b,maxbi)


    #methode de la bissectrice pour trouver les valeurs propres de T


    #methode de la puissance inverse pour trouver les vecteurs propres  de T
m = eigvals.shape[0]
Q = np.zeros((m,m))
for i in range(0,m):
        mu = eigvals[i]
        x = np.random.randint(10, size=(m)) #vecteur quelconque
        M = T - mu * np.eye(m)
        P, L, U = nla.lu (M)
        for k in range (0, 10) :
          v = np.dot(np.transpose(P),x) #P orthogonale (ici P != I)
          descente(L, v)      #ou v = nla.solve_triangular(L,v), lower=True)
          remonte(U, v)       #ou v = nla.solve_triangular(U, v), lower=False)
          v = (1 / nla.norm(v,2)) * v
          Q[:,i] = v




    Lambda = eigvals [n:2*n] #extraction des valeurs propres positives qui sont les valeurs singulières de A
    Q = Q [:,n:2*n]
    Q = eivalgs[:,n:2*n]
m,n = A.shape
Y = np.sqrt(2) * np.dot (P, Q)
newVt = np.transpose(Y[0:n,:])
newU = np.zeros ([m,n], dtype=np.float64)
newU[0:n,:] = Y[n:2*n,:]
newSigma = np.array (np.diag (Lambda), dtype=np.float64)    #B = np.dot(newU,np.dot(newSigma,newVt))
for i in range (n-1, -1, -1) :        #E4
  Q = reflecteur (m, VL[i])
  newU = np.dot (Q, newU)
for i in range (n-3, -1, -1) :
  Q = reflecteur (n, VR[i])
  newVt = np.dot (newVt,Q) #newA = np.dot (newU, np.dot (newSigma, newVt))
return newU, newSigma, newVt
                   

U, sigma, Vt = calcul_svd(A)
print("U =",U, "\nsigma =",sigma, "\nVt =",Vt)


