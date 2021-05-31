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
U, val_sing, Vt = nla.svd(A)
Sigma = np.diag(val_sing)


#Q5
v1 = Vt[0,:]
v2 = Vt[1,:]

P=np.dot(A, np.transpose(np.array([v1, v2])))

plt.scatter (P[:,0], P[:,1])
for i in range (0,m) :
    plt.annotate (lignes[i], P[i,:])

plt.scatter (0, 0, color='white')
plt.annotate ('', np.zeros(2), ha='center')

plt.show ()




