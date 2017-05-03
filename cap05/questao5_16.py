import numpy as np
from numpy import linalg as LA
from math import sqrt

def newton(f, jac, x0):
	
	b=-f(x0)
	print "det = ", LA.det(jac(x0))

	while LA.norm(b) > 1e-12:
		print "r = ", LA.norm(b)

		j = jac(x0)

		x = LA.solve(j, b)

		x0 = x0 + x
		b=-f(x0)
		
	print x0
	print LA.norm(b)


def f1(x):
	r = np.zeros(3)

	r[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 5

	r[1] = x[0] + x[1] - 1

	r[2] = x[0] + x[2] - 3

	return r


def jac1(x):

	j = np.zeros((3,3))

	j[0,0] = 2.0*x[0]
	j[0,1] = 2.0*x[1]
	j[0,2] = 2.0*x[2]

	j[1,0] = 1.0
	j[1,1] = 1.0
	j[1,2] = 0.0

	j[2,0] = 1.0
	j[2,1] = 0.0
	j[2,2] = 1.0

	return j


def f2(xv):
	r = np.zeros(3)
	x,y,z = xv[0],xv[1],xv[2]
	r[0] = 16*(x*x*x*x + y*y*y*y) + z*z*z*z - 16
	r[1] = x*x + y*y + z*z - 3
	r[2] = x*x*x - y

	return r


def jac2(xv):
	jac = np.zeros((3,3))
	x,y,z = xv[0],xv[1],xv[2]
	jac[0,0] = 16*4*x*x*x
	jac[0,1] = 16*4*y*y*y
	jac[0,2] = 4*z*z*z

	jac[1,0] = 2*x
	jac[1,1] = 2*y
	jac[1,2] = 2*z	

	jac[2,0] = 3.0*x*x
	jac[2,1] = -1.0
	jac[2,2] = 0.0

	return jac

x0 = np.array([1.0,1.0,1.0])
x0 = np.array([(1+sqrt(3))/2, (1-sqrt(3))/2, sqrt(3)])
newton(f1, jac1, x0)