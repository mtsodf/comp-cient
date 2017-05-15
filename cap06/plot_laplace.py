import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import os
import sys





os.system("./build/Laplaciano %s" % sys.argv[1])

calc = []
sol  = []
f = open("./saida.txt")
state = 1
for line in f:

	if state == 1:
		n = int(line.split()[2])
		state = 2
		f.next()

		continue

	if state == 2:
		for i in range(n-1):
			line = f.next()
			calc.extend([float(x) for x in line.strip().split()])

		state = 3

		f.next();
		continue

	if state == 3:
		for i in range(n-1):
			line = f.next()
			sol.extend([float(x) for x in line.strip().split()])
		state = 4

x = [(1.0/n)*i for i in range(1, n)]
y = [(1.0/n)*i for i in range(1, n)]

calc = np.array(calc).reshape(n-1,n-1)
sol = np.array(sol).reshape(n-1,n-1)

   
    

fig = plt.figure(1, (8,12))
    
my_cmap = plt.cm.get_cmap('seismic')
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(2, 1),
                     axes_pad=0.4,
                     add_all=True,
                     label_mode="L",
                     cbar_location="right",
                     cbar_mode="each"
                     )


ax = grid[0]
im = ax.imshow(calc, interpolation="none", cmap='seismic', origin="lower")
grid.cbar_axes[0].colorbar(im)

ax = grid[1]
im = ax.imshow(sol, interpolation="none", cmap='seismic', origin="lower")
grid.cbar_axes[1].colorbar(im)

plt.show()