#PROGRAMA QUE GENERA LAS GRÁFICAS PARA EL RECONOCIMIENTO DE LETRAS
#INPUT: Fichero obtenido de "reconocimiento.cpp" con los solapamientos "reconocimiento.txt"
#OUTPUT: Gráfica del reconocimiento de patrones "reconocimiento.png"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sol = np.loadtxt("reconocimiento.txt") #ABRIMOS FICHERO DE DATOS

x = sol[:,0]
solA = sol[:,1]
solB = sol[:,2]
solC = sol[:,3]

indexA = np.where(sol[:,0] == 1)
solAm = np.average(abs(sol[indexA,1])) #SOLAPAMIENTO MEDIO PARA EL PATRÓN "A"
indexB = np.where(sol[:,0] == 2)
solBm = np.average(abs(sol[indexB,2])) #SOLAPAMIENTO MEDIO PARA EL PATRÓN "B"
indexC = np.where(sol[:,0] == 3)
solCm = np.average(abs(sol[indexC,3])) #SOLAPAMIENTO MEDIO PARA EL PATRÓN "C"

my_colors = ["#293462", "#F24C4C", "#EC9B3B"]
df = pd.DataFrame( {'lab':['A', 'B', 'C'], 'val':[solAm,solBm,solCm]} )
df.plot.bar(x='lab', y='val', rot=0, color=my_colors, legend=None)

#OBTENEMOS LA GRÁFICA
plt.title("Reconocimiento de patrones")
plt.xlabel("Patrón inicial deformado")
plt.ylabel("Solapamiento medio final ($\overline{m}^{\mu}$)")
plt.tight_layout()
plt.savefig("reconocimiento.png", dpi=300)
plt.show()