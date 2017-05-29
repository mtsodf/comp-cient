'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math

def plotar_tempo(ax, t, n):


    # Make data.
    X = np.arange(1.0/n, 1.0, 1.0/n)
    Y = np.arange(1.0/n, 1.0, 1.0/n)
    X, Y = np.meshgrid(X, Y)
    R = 20*X * (1 - X) * Y * (1 - Y);
    Z = R * math.exp(-t)

    ax.set_title('Analitica')
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0.0, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))



def plot_3d(Z, n, t):

    fig = plt.figure(figsize=plt.figaspect(2.))
    fig.suptitle('t = %f' % t)
    ax = fig.add_subplot(2, 1, 1, projection='3d')
    ax.set_title('Calculada')

    # Make data.
    X = np.arange(1.0/n, 1.0, 1.0/n)
    Y = np.arange(1.0/n, 1.0, 1.0/n)
    X, Y = np.meshgrid(X, Y)

    Z = np.reshape(Z, (n-1,n-1))

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0.0, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    ax = fig.add_subplot(2, 1, 2, projection='3d')
    plotar_tempo(ax, t, n)



    fig.savefig('./solucao/t_%4.2f.png'%t)   # save the figure to file

    plt.close(fig)


def read_saida():
    f = open("saida.txt")

    n = 0;
    unknows = (n-1)*(n-1)
    valores = {}
    t = None
    tempos = []
    cont = 0
    for line in f:
        if n == 0:
            n = int(line)
        else:
            if t is None:
                t = float(line)
                tempos.append(t)
                valores[t] = []
                cont = 0
            else:
                valores[t].extend([float(x) for x in line.split()])
                cont += 1
                if cont == n-1:
                    t = None


    print "qtd = ", len(valores[tempos[0]])
    print n

    for t in tempos[:20]:
        plot_3d(valores[t], n, t)
        #plotar_tempo(t,n)


read_saida()

#plotar_tempo(0.0)