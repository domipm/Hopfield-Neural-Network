#PROGRAMA QUE GENERA LAS GRÁFICAS PARA RECUPERACIÓN DE LA MEMORIA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from colour import Color

sol = np.loadtxt("memoria.txt") #ABRIMOS FICHERO DE DATOS

x = sol[:,0]
y = sol[:,1]

init = Color("#293462")
colors = list(init.range_to(Color("#F24C4C"),50))
colors = [color.rgb for color in colors]

df = pd.DataFrame({'lab':x, 'val':y}).astype(int)
df.plot.bar(x='lab', y='val', rot=0, legend=None, color=colors)

#OBTENEMOS LA GRÁFICA INICIAL
plt.title("Recuperación de la memoria")
plt.xlabel("Número de patrones almacenados $(n_{pat}$)")
plt.ylabel("Patrones recordados")
plt.xlim((-0.5,29.5))
plt.ylim((0,1.05))
plt.xticks(range(0,30,2))
plt.yticks(range(0,2,1))
plt.ticklabel_format(axis="y", style="plain")
plt.tight_layout()
plt.savefig("memoria.png", dpi=300)
plt.show()

plt.close()

yn = y/20**2

init = Color("#293462")
colors = list(init.range_to(Color("#F24C4C"),50))
colors = [color.rgb for color in colors]

df = pd.DataFrame({'lab':x.astype(int), 'val':yn})
df.plot.bar(x='lab', y='val', rot=0, legend=None, color=colors)

#OBTENEMOS LA GRÁFICA PARA LA FRACCIÓN
plt.title("Fracción máxima patrones recordados")
plt.xlabel("Número de patrones almacenados $(n_{pat}$)")
plt.ylabel("Fracción $(\\alpha_c)$")
plt.xlim((-0.5,29.5))
plt.ylim((0,0.005))
plt.xticks(range(0,30,2))
plt.ticklabel_format(axis="y", style="plain")
plt.tight_layout()
plt.savefig("fraccion.png", dpi=300)
plt.show()