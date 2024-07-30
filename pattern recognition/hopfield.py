import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as col

plt.rcParams["animation.convert_path"] = "/usr/bin/convert"

fname = "output"

data = np.loadtxt(fname+".txt")

color_dict = {0: '#2D4263', 1: '#C84B31'}
label_dict = {0: 'off', 1: 'on'}
imax = max(label_dict)
imin = min(label_dict)
cmap = col.ListedColormap(color_dict.values())

dim = int(len(data[0]))
iter = int(len(data[:,0])/(dim))

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

plt.title("Pattern Recognition")
plt.colorbar(mat,ticks=(0,1))
#plt.show()
ani.save(fname+".gif", writer=animation.PillowWriter(fps=3), dpi=300)

plt.close()

exit()

sol = np.loadtxt("solapamiento.txt")
npat = len(sol[0,:])
x = range(0,iter)

plt.title("Solapamiento patrones")
plt.xlabel("Pasos Monte Carlo")
plt.ylabel("Solapamiento")
plt.xticks(range(0,22,2))
plt.ylim((-1.05,1.05))
plt.tight_layout()

for i in range(0, npat):
    plt.plot(x, sol[:,i], label="Patrón " + str(i+1))

plt.legend()
plt.savefig("solapamiento.png", dpi=300)
plt.show()