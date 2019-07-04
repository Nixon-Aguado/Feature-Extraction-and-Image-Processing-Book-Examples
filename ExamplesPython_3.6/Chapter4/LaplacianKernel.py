'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
LaplacianKernel: Creates a Laplacian kernel of arbitrary size
'''

# Set module functions
from ImageUtilities import createImageF
from PrintUtilities import printImageRangeF
from ImagePropertiesUtilities import imageMaxMin
from PlotUtilities import plotSurface

# Math and iteration
from math import exp
from timeit import itertools

'''
Parameters:
    kernelSize = Size of the kernel
    sigma = Standard deviation of the kernel
'''
kernelSize = 15
sigma = 1.5

# To store kernel
kernelLaplacian = createImageF(kernelSize, kernelSize)

# Create kernel
s2Inv = 1.0 / (sigma * sigma)
kernelCentre = (kernelSize - 1) / 2

# Generate kernel values
sumValues = 0.0
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):        
    nx2 = float(x-kernelCentre) * float(x-kernelCentre)
    ny2 = float(y-kernelCentre) * float(y-kernelCentre)
    s = 0.5 * (nx2 + ny2) * s2Inv
    
    kernelLaplacian[y,x] = - s2Inv * s2Inv * (1.0 - s) * exp(-s)
    sumValues += kernelLaplacian[y,x]

# Normalize
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    kernelLaplacian[y,x] /= sumValues
                                   
# Print kernel
printImageRangeF(kernelLaplacian, [0, kernelSize-1], [0, kernelSize-1], ' 8.2f')

# Plot surface
maxValue, minValue = imageMaxMin(kernelLaplacian)
plotSurface(kernelLaplacian, [minValue, maxValue], 1)
