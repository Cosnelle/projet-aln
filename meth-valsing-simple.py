import numpy as np
import math
import scipy.linalg as nla
import matplotlib.pyplot as plt

# Affichage plus agréable
np.set_printoptions(linewidth=240)

#A=np.array([[11.,32.,3.],[24.,5.,56.],[27.,8.,9.],[10.,11.,12.]])
A = np.array ([
    [3, 97.8, 119, 2.6, 7.4, 0.2, 1.8, 2.2, 53, 34, 28.3, 21],
    [0.2, 80.5, 121, 2.1, 5.3, 0.5, 1.4, 1.4, 47.4, 20.5, 27.8, 19.9],
    [0.4, 6.1, 67, 4.2, 9.9, 3, 1.5, 2.5, 10.6, 21, 15.1, 23.1],
    [7.3, 106.4, 129, 1.9, 14.7, 1.1, 2.1, 3.7, 45.9, 12.5, 22.1, 29.9],
    [4.6, 170.6, 79, 1, 26.4, -4.4, 1.4, 11.8, 27.6, 23, 26, 31],
    [6.7, 69.3, 98, 2.4, 26.2, -1.4, 1.4, 4.9, 32.7, 35, 22.7, 27],
    [3.7, 86, 108, 2.2, 10.6, 0.1, 2, 2, 45.4, 34, 30.8, 19.3],
    [2.1, 120.7, 100, 3.3, 11.7, -1, 1.4, 4.6, 35.6, 33, 27.8, 24.5],
    [4.5, 71.1, 94, 3.1, 14.7, -3.5, 1.4, 12, 20.7, 10, 18.4, 23.5],
    [0.9, 18.3, 271, 2.9, 5.3, 0.5, 1.5, 2.3, 50.7, 22, 20.1, 16.8],
    [2.9, 70.9, 85, 3.2, 7, 1.5, 1.5, 4, 25, 35, 18.9, 21.4],
    [3.6, 65.5, 131, 2.8, 6, -0.6, 1.8, 1.6, 55.7, 34, 28.4, 15.7],
    [2.5, 72.4, 129, 2.6, 4.9, 0.7, 1.4, 1.7, 45.6, 34, 28.2, 16.9],
    [4.9, 108.1, 77, 2.8, 17.6, -1.9, 1.4, 6, 14.8, 25, 24.3, 24.4],
    [5.1, 46.9, 84, 2.8, 10.2, -2, 1.6, 4.9, 18.5, 20, 21.5, 19.3],
    [3.3, 43.3, 73, 3.7, 14.9, 1.1, 1.5, 2.1, 10.1, 19, 16, 20.6],
    [1.5, 49, 114, 3.2, 7.9, 0.3, 1.8, 1.6, 46.8, 26, 26.3, 17.9]],
    dtype=np.float64)

lignes = ['BE', 'DE', 'EE', 'IE', 'EL', 'ES', 'FR', 'IT', 'CY', 'LU', 'MT', 'NL', 'AT', 'PT', 'SL', 'SK', 'FI']

m,n = A.shape



def centrer(A):
  m,n = A.shape
  for j in range (0,n):
    x = A[:,j]
    Moy = np.mean(x)
    for i in range (0,m):
       A[i][j] = A[i][j]-Moy
  return A

centrer(A)


def reduit(A):
    m,n = A.shape
    for j in range (0,n):
        x = A[:,j]
        E = np.std(x)
        for i in range (0,m):
            A[i][j] = A[i][j]/E
    return A

reduit(A)


#Q4
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

   
    #methode de la bissectrice pour trouver les valeurs propres de T


    #methode de la puissance inverse pour trouver les vecteurs propres  de T
    m = eigvals.shape[0]
    Q = np.zeros((m,m))
    for i in range(0,m):
        mu = eigvals[i]
        x = np.random.randint(10, size=(m)) #vecteur quelconque
        M = T - mu * np.eye(m)
        Pb, L, U = nla.lu (M)
        for k in range (0, 10) :
          v = np.dot(np.transpose(Pb),x) #Pb orthogonale (ici Pb != I)
          descente(L, v)      #ou v = nla.solve_triangular(L,v), lower=True)
          remonte(U, v)       #ou v = nla.solve_triangular(U, v), lower=False)
          v = (1 / nla.norm(v,2)) * v
          Q[:,i] = v


    Lambda = eigvals [n:2*n]
    Q = Q [:,n:2*n]
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
      newVt = np.dot (newVt, Q)     #newA = np.dot (newU, np.dot (newSigma, newVt))
    return newU, newSigma, newVt


U, sigma, Vt = calcul_svd(A)


#Q5
v1 = Vt[n-1,:]
v2 = Vt[n-2,:]

P=np.dot(A, np.transpose(np.array([v1, v2])))

plt.scatter (P[:,0], P[:,1])
for i in range (0,m) :
    plt.annotate (lignes[i], P[i,:])

plt.scatter (0, 0, color='white')
plt.annotate ('', np.zeros(2), ha='center')

plt.show ()




