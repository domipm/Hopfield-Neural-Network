# Prepare any image to use as input matrix for the model

import os
from PIL import Image, ImageEnhance
import numpy as np

contrast = 8
newsize = 128 # Number of pixels (per side, square image)
fname = "moon"

img = Image.open(fname+".jpg")                                  # Open image
enhancer = ImageEnhance.Contrast(img)
img = enhancer.enhance(contrast)                                # Contrast adjustment
img = img.resize((newsize, newsize),resample=Image.BILINEAR)
img = img.convert('1')                                          # Black and white
img = img.resize((newsize,newsize),Image.NEAREST)               # Restore original size

img.save(fname+"_result.png")   # Save new image

data = open(fname+".txt", "w")  # Create output file

pixel = np.zeros((newsize,newsize), dtype=int) # Binary pixel array of photo (0,1)

# Obtain binary value for each pixel
pix = img.load()
for x in range(newsize):
    for y in range(newsize):
        if (pix[x,y] == 255): pixel[x,y] = 1
        else: pixel[x,y] = 0

np.savetxt(data, np.transpose(pixel), fmt="%i") # Save pixel array on new file

data.close() # Close file