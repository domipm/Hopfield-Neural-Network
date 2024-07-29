import os
from PIL import Image, ImageEnhance
import numpy as np

contrast = 8
newsize = 128 #NUMERO DE PIXELES 
fname = "moon"

img = Image.open(fname+".jpg") #ABRIMOS FOTO
enhancer = ImageEnhance.Contrast(img)
img = enhancer.enhance(contrast) #CONTRASTE
img = img.resize((newsize, newsize),resample=Image.BILINEAR)
img = img.convert('1') #BLANCO Y NEGRO
img = img.resize((newsize,newsize),Image.NEAREST) #DEVOLVEMOS AL TAMAÑO ORIGINAL

img.save(fname+"_result.png") #GUARDAMOS FOTO NUEVA

data = open(fname+".txt", "w") #CREAMOS FICHERO DE SALIDA

pixel = np.zeros((newsize,newsize), dtype=int) #ARRAY DE PIXELES DE LA FOTO (0, 1)

#OBTENEMOS VALOR BINARIO DE CADA PIXEL
pix = img.load()
for x in range(newsize):
    for y in range(newsize):
        if (pix[x,y] == 255): pixel[x,y] = 1
        else: pixel[x,y] = 0

np.savetxt(data, np.transpose(pixel), fmt="%i") #GUARDAMOS ARRAY DE PIXELES EN EL ARCHIVO

data.close() #CERRAMOS FICHERO