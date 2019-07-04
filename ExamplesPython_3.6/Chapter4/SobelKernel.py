'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
SobelKernel: Generate a Sobel kernel of arbitrary size 
'''

# Set module functions
from ImageUtilities import createImageUV
from PrintUtilities import printImageRangeF

# Math and iteration
from math import factorial
from timeit import itertools 

'''
Parameters:
    kernelSize = Size of the kernel
'''
kernelSize = 5

# Pascal kernels pascal2 is a shifted version of pascal1
pascal1 = createImageUV(kernelSize, kernelSize)
pascal2 = createImageUV(kernelSize, kernelSize)

smooth = createImageUV(kernelSize, kernelSize)
sobel = createImageUV(kernelSize, kernelSize)

# Create kernel
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
        
    # Smooth
    smooth[y,x,0] = factorial(kernelSize - 1) /                           \
                    (factorial(kernelSize - 1 - x) * factorial(x))
    smooth[y,x,1] = factorial(kernelSize - 1) /                           \
                    (factorial(kernelSize - 1 - y) * factorial(y))
    
    # Pascal 
    if (kernelSize - 2 - x >= 0):
        pascal1[y,x,0] = factorial(kernelSize - 2) /                      \
                         (factorial(kernelSize - 2 - x) * factorial(x))
    
    if (kernelSize - 2 - y >= 0):    
        pascal1[y,x,1] = factorial(kernelSize - 2) /                      \
                         (factorial(kernelSize - 2 - y) * factorial(y))
        
    # Pascal shift to the right
    xp = x - 1
    if (kernelSize - 2 - xp >= 0 and xp >= 0):
        pascal2[y,x,0] = factorial(kernelSize - 2) /                      \
                         (factorial(kernelSize - 2 - xp) * factorial(xp))
        
    yp = y - 1
    if (kernelSize - 2 - yp >= 0 and yp >= 0):
        pascal2[y,x,1] = factorial(kernelSize - 2) /                      \
                         (factorial(kernelSize - 2 - yp) * factorial(yp))
        
    # Sobel    
    sobel[y,x,0] = smooth[y,x,1] * (pascal1[y,x,0] - pascal2[y,x,0])
    sobel[y,x,1] = smooth[y,x,0] * (pascal1[y,x,1] - pascal2[y,x,1])                                       
        
# Print pixel's values of the kernel
printImageRangeF(smooth[:,:,0], [0, kernelSize-1], [0, kernelSize-1], '2.0f')
printImageRangeF(smooth[:,:,1], [0, kernelSize-1], [0, kernelSize-1], '2.0f')

printImageRangeF(pascal1[:,:,0], [0, kernelSize-1], [0, kernelSize-1], '2.0f')
printImageRangeF(pascal1[:,:,1], [0, kernelSize-1], [0, kernelSize-1], '2.0f')

printImageRangeF(pascal2[:,:,0], [0, kernelSize-1], [0, kernelSize-1], '2.0f')
printImageRangeF(pascal2[:,:,1], [0, kernelSize-1], [0, kernelSize-1], '2.0f')

printImageRangeF(sobel[:,:,0], [0, kernelSize-1], [0, kernelSize-1], '4.0f')
printImageRangeF(sobel[:,:,1], [0, kernelSize-1], [0, kernelSize-1], '4.0f')
