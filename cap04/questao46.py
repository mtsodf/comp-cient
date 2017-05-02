import numpy as np
from numpy import linalg as LA

def Lanczos(n):

	B = np.random.rand(n,n)

	Q, R = np.linalg.qr(B)

	D = np.zeros((n,n))

	for i in range(n):
		D[i,i] = i + 1


	A = np.dot(D, np.transpose(Q))
	A = np.dot(Q, A)


	w, q = LA.eig(A)

	print "Autovalores ", np.sort(w)

	r = np.random.rand(n)
	r = r/LA.norm(r)

	beta = np.zeros(n+1)
	alpha = np.zeros(n+1)

	beta[0] = LA.norm(r)

	qant = np.zeros(n)

	eigens = []
	for i in range(1,n+1):
		q = r / beta[i-1]
		u = np.dot(A, q)
		r = u - beta[i-1] * qant
		alpha[i] = np.dot(q, r)
		r = r - alpha[i]*q
		beta[i] = LA.norm(r)

		qant = q

		T = np.zeros((i, i))

		for j in range(i):
			T[j,j] = alpha[j+1]

		for j in range(i-1):
			T[j, j+1] = beta[j+1]
			T[j+1, j] = beta[j+1]

		w, v = LA.eig(T)

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
	ax.set_title("n = %d" % n)
	plt.savefig('n_%d.png'%n)


Lanczos(10)
Lanczos(20)
Lanczos(30)
Lanczos(40)
Lanczos(50)