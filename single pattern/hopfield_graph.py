#PROGRAMA QUE GENERA EL GIF DE LA EVOLUCIÓN DE LA MATRIZ SISTEMA Y LA GRÁFICA DE SOLAPAMIENTO
#INPUT: Fichero entrada "output.txt" obtenido de "hopfield.cpp"
#OUTPUT: Animación "output.gif"
#OUTPUT: Gráfica de solapamiento en función de los pasos Monte Carlo "solapamiento.png"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as col

plt.rcParams["animation.convert_path"] = "/usr/bin/convert"

fname = "output"

data = np.loadtxt(fname+".txt") #FICHERO DE MATRICES
print(data)

color_dict = {0: '#2D4263', 1: '#C84B31'}
label_dict = {0: 'off', 1: 'on'}
imax = max(label_dict)
imin = min(label_dict)
cmap = col.ListedColormap(color_dict.values())

dim = int(len(data[0])) #DIMENSION DE LA MATRIZ
iter = int(len(data[:,0])/(dim)) #NUMERO DE ITERACIONES

#MATRIZ
grid = np.zeros((dim,dim))

def update(k):
    global grid
    newGrid = grid.copy()
    for i in range(dim):
        for j in range(dim):
            newGrid[i][j] = data[i+k*dim][j]
    mat.set_data(grid)
    grid = newGrid
    return [mat]

fig, ax = plt.subplots()
mat = ax.matshow(grid,cmap=cmap,interpolation='nearest',vmin=imin,vmax=imax)
ani = animation.FuncAnimation(fig, update, frames=np.arange(0, iter, 1), interval=1, repeat=False)

plt.title("Hopfield Neural Network")
plt.tight_layout()
plt.colorbar(mat,ticks=(0,1))

plt.show()

ani.save(fname+".gif", writer=animation.PillowWriter(fps=3), dpi=300)

plt.close()

exit()

sol = np.loadtxt("solapamiento.txt")
x = range(0,iter)

xzoom = np.arange(0,5,1)
yzoom = sol[0:5]

plt.title("Solapamiento según pasos Monte Carlo")
plt.xlabel("Pasos Monte Carlo")
plt.ylabel("Solapamiento")
plt.ylim((-1.05,1.05))
plt.tight_layout()
plt.grid(False)
plt.xticks(range(0,21,2))
plt.plot(x, sol, label="Patrón 1")

axes = plt.axes([.4, .2, .5, .5])
axes.plot(xzoom, yzoom, c="green")
axes.grid(False)

plt.savefig("solapamiento.png", dpi=300)
plt.show()