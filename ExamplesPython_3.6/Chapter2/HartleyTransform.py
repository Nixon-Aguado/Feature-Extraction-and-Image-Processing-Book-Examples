'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 2
HartleyTransform: Compute the Hartley transform of an image
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from ImageOperatorsUtilities import imageLogF
from PrintUtilities import printProgress

# Math functions and iteration
from math import sin, cos, pi
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter2/Input/"
imageName = "Giraffe.png" 

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create coefficients Image. Maximum frequency according to sampling
maxFreqW = int(width  / 2)
maxFreqH = int(height / 2)
numCoeffW = 1 + 2 * maxFreqW
numCoeffH = 1 + 2 * maxFreqH
coeff = createImageF(numCoeffW ,numCoeffH)

# Adjust the size of the data to be even
m = float(width)
n = float(height)
if width % 2 == 0:
    m = width + 1.0
if height % 2 == 0:
    n = height + 1.0  
    
# Fundamental frequency
ww = (2.0 * pi) / m
wh = (2.0 * pi) / n

# Compute values
for u in range(-maxFreqW, maxFreqW + 1):
    printProgress(u + maxFreqW, numCoeffW)
    entryW = u + maxFreqW
    for v in range(-maxFreqH, maxFreqH + 1):
        entryH = v + maxFreqH 
        for x,y in itertools.product(range(0, width), range(0, height)):
            coeff[entryH, entryW] += inputImage[y,x] *                                      \
                 (cos(x * ww * u) + sin(x * ww * u)) * (cos(y * wh * v) + sin(y * wh * v)) 

# Include scale
for u in range(-maxFreqW, maxFreqW + 1):
    printProgress(u + maxFreqW, numCoeffW)
    entryW = u + maxFreqW
    for v in range(-maxFreqH, maxFreqH + 1):
        entryH = v + maxFreqH 
        coeff[entryH, entryW] /= m*n 

# Show transform in log form. The function converts negative values to positive, 
# so it is similar to the power 
coeffLog = imageLogF(coeff)
showImageF(coeffLog)
