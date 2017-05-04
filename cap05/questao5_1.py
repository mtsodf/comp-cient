import numpy as np
from numpy import linalg as LA
from math import sqrt
import matplotlib.pyplot as plt

def arrow(ax, x1, y1, x2, y2):
	hl = 0.1
	dx = x2-x1
	dy = y2-y1

	if(x2 > x1):
		dx -= hl
	if(x1 > x2):
		dx += hl

	if(dy > 0):
		dy -= hl

	if(dy<0):
		dy += hl

	ax.arrow(x1, y1, dx, dy, head_width=0.05, head_length=hl, fc='k', ec='k')

def fixed_point(filename, f, x0, tol=1e-5, xmin=-2, xmax=2):

	iters = 0
	xs = []
	while (iters<10 and abs(f(x0))>tol):
		print iters
		xs.append(x0)
		x0 = f(x0)
		print x0
		iters += 1

	if(abs(f(x0))<tol):
		print "Convergiu para %f" % x0

	else:
		"Nao convergiu"


	
	fig, ax = plt.subplots()
	ax.arrow(0, xmin, 0, xmax-xmin, head_width=0.05, head_length=0.1, fc='k', ec='k')
	ax.arrow(xmin, 0, xmax-xmin, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

	linex = np.linspace(xmin,xmax,num=30)
	ax.plot(linex,linex, "b-")
	ax.plot(linex,[f(x) for x in linex], "r-")
	ax.plot(xs[0], 0, "ro")
	ax.plot(xs, [f(x) for x in xs], "ro")
	ax.grid(True)



	arrow(ax, xs[0],0,xs[0],f(xs[0]))

	for i in range(1, len(xs) - 1):
		arrow(ax, xs[i-1],f(xs[i-1]), xs[i], xs[i])
		arrow(ax, xs[i], xs[i], xs[i],f(xs[i]))


		

	plt.savefig(filename)


fa = lambda x: x*x-2
fb = lambda x: sqrt(x+2)
fc = lambda x: 1 + 2/x
fd = lambda x: (x*x + 2)/(2*x-1)

fixed_point("item_a", fa, -0.5)
fixed_point("item_b", fb, -0.5, xmin=-1, xmax=3)
fixed_point("item_b_2", fb,  10.0, xmin=0.0, xmax=12.0)
fixed_point("item_c", fc, 5.0, xmin=0.0, xmax=6.0)
fixed_point("item_d", fd, 10.0, xmin=-0.0, xmax=10.0)
fixed_point("item_d_2", fd, -4.0, xmin=-10.0, xmax=0.0)
