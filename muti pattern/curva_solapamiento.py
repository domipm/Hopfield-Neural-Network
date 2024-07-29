#PROGRAMA PARA REALIZAR LA GRÁFICA DE LA CURVA DE SOLAPAMIENTO FRENTE A TEMPERATURA

import numpy as np
import scipy.stats as stats
from scipy import optimize
import matplotlib.pyplot as plt

sol = np.loadtxt("curva_solapamiento.txt") #ABRIMOS FICHERO DE DATOS

#OBTENEMOS LA GRÁFICA
plt.title("Curva de solapamiento frente a temperatura")
plt.xlabel("Temperatura ($T$)")
plt.ylabel("Solapamiento ($m^{\mu}$)")
plt.ylim((-0.2,1.05))
plt.tight_layout()
plt.grid(False)
plt.xscale("log")

plt.plot(sol[:,0], sol[:,1], marker=".", linestyle="-", alpha=1, markersize=5, label="$m^1$")
plt.plot(sol[:,0], sol[:,2], marker=".", linestyle="-", alpha=1, markersize=5, label="$m^2$")
plt.plot(sol[:,0], sol[:,3], marker=".", linestyle="-", alpha=1, markersize=5, label="$m^3$")
plt.legend()
plt.savefig("curva_solapamiento.png", dpi=300)
plt.show()