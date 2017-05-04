import numpy as np
from numpy import linalg as LA
from math import sqrt, exp

def newton(f, jac, x0):
	
	b=-f(x0)
	

	print "%12s    |%12s    |%12s" % ("iteracao", "|r|", "det(jac)")
	iter = 0
	while LA.norm(b) > 1e-12 and iter<1000:

		j = jac(x0)
		print "%12d    |%12.4e    |%12.4e" % (iter, LA.norm(b), LA.det(j))

		x = LA.solve(j, b)

		x0 = x0 + x
		b=-f(x0)

		iter += 1

	print "%12d    |%12.4e    |" % (iter, LA.norm(b))

	print "x = ", x0

#ITEM A
def f1(x):
	r = np.zeros(2)
	x1 = x[0]
	x2 = x[1]
	
	r[0] = x1 + 5*x2*x2 - x2*x2*x2 - 2*x2 - 13
	r[1] = x1 + x2*x2 + x2*x2*x2 - 14*x2 - 29

	return r

def jac1(x):
	x1 = x[0]
	x2 = x[1]

	j = np.zeros((2,2))

	j[(0,0)] = 1.0
	j[(1,0)] = 1.0

	j[(0,1)] = 10*x2 - 3*x2*x2 -2
	j[(1,1)] = 2*x2 + 3*x2*x2 - 14

	return j 



# ITEM B
def f2(x):
	r = np.zeros(3)

	r[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 5

	r[1] = x[0] + x[1] - 1

	r[2] = x[0] + x[2] - 3

	return r


def jac2(x):

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

#ITEM C

def f3(x):
	x1,x2,x3,x4 = x

	r = np.zeros(4)

	r[0] = x1 + 10*x2
	r[1] = sqrt(5)*(x3 - x4)
	r[2] = (x2 - x3) * (x2 - x3)
	r[3] = sqrt(10)*(x1-x4)*(x1-x4)

	return r

def jac3(x):
	x1,x2,x3,x4 = x
	j = np.zeros((4,4))

	j[0,0] = 1.0
	j[0,1] = 10.0
	j[0,2] = 0.0
	j[0,3] = 0.0

	j[1,0] = 0.0
	j[1,1] = 0.0
	j[1,2] = sqrt(5)
	j[1,3] = sqrt(5)

	j[2,0] = 0.0
	j[2,1] = 2*(x2-x3)
	j[2,2] = -2*(x2-x3)
	j[2,3] = 0.0	

	j[3,0] = sqrt(10)*2*(x1-x4)
	j[3,1] = 0.0
	j[3,2] = 0.0
	j[3,3] = -sqrt(10)*2*(x1-x4)	

	return j

#ITEM D
def f4(x):
	x1,x2 = x

	r = np.zeros(2)

	r[0] = x1
	r[1] = 10*x1/(x1+0.1) + 2*x2*x2

	return r

def jac4(x):
	x1,x2 = x

	j = np.zeros((2,2))

	j[0,0] = 1
	j[0,1] = 0

	j[1,0] = pow((x1+0.1), -2)
	j[1,1] = 4*x2

	return j

def f5(x):
	x1, x2 = x

	r = np.zeros(2)

	r[0] = 10000*x1*x2 - 1

	r[1] = exp(-x1) + exp(-x2) - 1.0001

	return r

def jac5(x):
	j = np.zeros((2,2))

	x1,x2 = x

	j[0,0] = 10000*x2
	j[0,1] = 10000*x1

	j[1,0] = -exp(-x1)
	j[1,1] = -exp(-x2)

	return j



def fextra(xv):
	r = np.zeros(3)
	x,y,z = xv[0],xv[1],xv[2]
	r[0] = 16*(x*x*x*x + y*y*y*y) + z*z*z*z - 16
	r[1] = x*x + y*y + z*z - 3
	r[2] = x*x*x - y

	return r


def jacextra(xv):
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



print "ITEM A"
x0 = np.array([15.0, -2.0])
newton(f1, jac1, x0)


print "ITEM B"
x0 = np.array([(1+sqrt(3))/2, (1-sqrt(3))/2, sqrt(3)])
newton(f2, jac2, x0)


print "ITEM C"
try:
	x0 = np.array([1.0, 2.0, 1.0, 1.00])
	newton(f3, jac3, x0)
except Exception as e:
	print "ITEM C NAO CONVERGIU"

print "ITEM D"
try:
	x0 = np.array([1.8, 0.00])
	newton(f4, jac4, x0)
except Exception as e:
	print "ITEM D NAO CONVERGIU"

print "ITEM E"
try:
	x0 = np.array([0.0, 1.0])
	newton(f5, jac5, x0)
except Exception as e:
	print "ITEM E NAO CONVERGIU"