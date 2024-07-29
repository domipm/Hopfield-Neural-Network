#PROGRAMA PARA REALIZAR LA GRÁFICA DE LA CURVA DE SOLAPAMIENTO FRENTE A TEMPERATURA
#INPUT: Fichero "solapamientot.txt" obtenido del programa "hopfield_sol.cpp"
#OUTPUT: Gráfica "solapamientot.png" de la curva de solapamiento

import numpy as np
import scipy.stats as stats
from scipy import optimize
import matplotlib.pyplot as plt

#FUNCIÓN DE AJUSTE
def f(x, A, B, C, D):
    return A*np.tanh(B*x+C) + D

sol = np.loadtxt("solapamientot.txt") #ABRIMOS FICHERO DE DATOS

popt, pcov = optimize.curve_fit(f,sol[:,0],sol[:,1]) #REALIZAMOS AJUSTE

#ESCRIBIMOS EN PANTALLA LOS PARÁMETROS DE AJUSTE
print(popt)
print(pcov)

#CALCULAMOS Y MOSTRAMOS EN PANTALLA EL PARÁMETRO CHI^2
chisq = np.sum((f(sol[:,0], *popt)-sol[:,1])**2)
chisq = chisq
print(chisq)

#OBTENEMOS LA GRÁFICA
plt.title("Curva de solapamiento frente a temperatura")
plt.xlabel("Temperatura ($T$)")
plt.ylabel("Solapamiento ($m^1$)")
plt.ylim((-0.05,1.05))
plt.tight_layout()
plt.grid(False)
plt.xscale("log")
plt.vlines(0.027996, ymin=-0.05, ymax=1.05, color="orange", alpha=0.5, label="Punto crítico ($T_c$)")
plt.plot(sol[:,0], sol[:,1], marker=".", color="red", linestyle="", label="Solapamiento ($m^1$)")
plt.plot(sol[:,0], sol[:,1], marker="", color="red", linestyle="-", alpha=0.5)
plt.plot(sol[:,0], f(sol[:,0], *popt), color="black", label="Ajuste $f(T)$", linestyle="--")
plt.legend()
plt.savefig("solapamientot.png", dpi=300)
plt.show()