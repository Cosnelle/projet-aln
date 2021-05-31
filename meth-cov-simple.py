import numpy as np 
import math
import scipy.linalg as nla
import matplotlib.pyplot as plt

# Affichage plus agréable
np.set_printoptions(linewidth=240)

#A=np.array([[11.,32.,3.],[24.,5.,56.],[27.,8.,9.],[10.,11.,12.],[1,2,4]])
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

#Q1
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


#Q3
C = 1/(n-1)*np.dot(np.transpose(A),A)


#4
 

#Algorithme QR pour trouver une approximation des valeurs propres 
T = C
for k in range (0,100) :
    Q, R = nla.qr (T)
    T = np.dot (R, Q)

#sur la diagonale de T on a une approximation des valeurs propres de C




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

#Methode de la puissance inverse pour trouver les vecteurs propres 
#version reel classique

m = C.shape[0]
Q = np.zeros((m,m))
for i in range(0,m):
  mu = T[i,i] 
  x = np.random.randint(10, size=(m)) #vecteur quelconque
  M = C - mu * np.eye(m)
  P, L, U = nla.lu (M)
  for k in range (0, 10) :
       v = np.dot(np.transpose(P),x) #P orthogonale (ici P != I)
       descente(L, v)      #ou v = nla.solve_triangular(L,v), lower=True)
       remonte(U, v)       #ou v = nla.solve_triangular(U, v), lower=False)
       v = (1 / nla.norm(v,2)) * v
       Q[:,i] = v

#la matrice Q contient les vecteurs propres de C   

#verification
#C = np.dot(Q,np.dot(T,np.transpose(Q)))


#5
qi = Q[:,0]
qj = Q[:,1]


#6
P=np.dot(A, np.transpose(np.array([qi, qj])))


#7
m = A.shape[0]
plt.scatter (P[:,0], P[:,1])
for i in range (0,m) :
    plt.annotate (lignes[i], P[i,:])

plt.scatter (0, 0, color='white')
plt.annotate ('', np.zeros(2), ha='center')

plt.show ()




