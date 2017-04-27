import numpy as np
from numpy import linalg as LA



n = 100

B = np.random.rand(n,n)

Q, R = np.linalg.qr(B)

D = np.zeros((n,n))

for i in range(n):
	D[i,i] = i + 1


A = np.dot(D, np.transpose(Q))
A = np.dot(Q, A)


w, q = LA.eig(A)

print "Autovalores ", np.sort(w)

r = []
r.append(np.random.rand(n))
r[0] = r[0]/LA.norm(r)

beta = np.zeros(n+1)
alpha = np.zeros(n+1)

beta[0] = LA.norm(r)

qant = np.zeros(n)

eigens = []
for i in range(1,n+1):
	q = r[i-1] / beta[i-1]
	u = np.dot(A, q)
	r.append(u - beta[i-1] * qant)
	alpha[i] = np.dot(q, r[i])
	r[i] = r[i] - alpha[i]*q
	beta[i] = LA.norm(r[i])

	qant = q

	T = np.zeros((i, i))

	for j in range(i):
		T[j,j] = alpha[j+1]

	for j in range(i-1):
		T[j, j+1] = beta[j+1]
		T[j+1, j] = beta[j+1]

	w, q = LA.eig(T)

	eigens.append(w)



print "Autovalores ", np.sort(w)

x = []
y = []

for i in range(len(eigens)):
	for j in range(len(eigens[i])):
		y.append(i)
		x.append(eigens[i][j])

import matplotlib.pyplot as plt


fig, ax = plt.subplots()
ax.plot(x, y, 'ro')

plt.show()