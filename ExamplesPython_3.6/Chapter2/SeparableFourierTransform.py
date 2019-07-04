'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 2
SeparableFourierTransform: Compute the Fourier transform of an image using the separable formulation
                           Display the magnitude and phase and reconstruct image
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from ImageOperatorsUtilities import imageLogF
from PrintUtilities import printProgress

# Iteration and Math functions
from math import sin, cos, pi, sqrt, atan2
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter2/Input/"
imageName = "Square.png"

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create coefficients Image. Two floats to represent a complex number
# Maximum frequency according to sampling
maxFreqW = int(width /2)
maxFreqH = int(height/2)
numCoeffW = 1 + 2 * maxFreqW
numCoeffH = 1 + 2 * maxFreqH
coeff = createImageF(numCoeffW ,numCoeffH , 2)

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

# Fourier Transform
for u in range(-maxFreqW, maxFreqW + 1):
    printProgress(u + maxFreqW, numCoeffW)
    entryW = u + maxFreqW 
    
    for v in range(-maxFreqH, maxFreqH + 1):
        entryH = v + maxFreqH
        coeff[entryH, entryW] = [0, 0]
        
        for x in range(0, width):
            sumY = [0, 0] 
            
            for y in range(0, height):
                sumY[0] += inputImage[y,x] * cos(y * wh * v)
                sumY[1] += inputImage[y,x] * sin(y * wh * v)
            coeff[entryH, entryW][0] += sumY[0] * cos(x * ww * u) - sumY[1] * sin(x * ww * u)
            coeff[entryH, entryW][1] -= cos(x * ww * u) * sumY[1] + sin(x * ww * u) * sumY[0]

for kw in range(-maxFreqW, maxFreqW + 1):
    printProgress(kw + maxFreqW, numCoeffW)
    entryW = kw + maxFreqW 
    for kh in range(-maxFreqH, maxFreqH + 1):
        entryH = kh + maxFreqH
        coeff[entryH, entryW][0] *= m*n 
        coeff[entryH, entryW][1] *= m*n

# Reconstruction
reconstruction = createImageF(width, height)
for u in range(-maxFreqW, maxFreqW + 1):
    printProgress(u + maxFreqW, numCoeffW)
    entryW = u + maxFreqW 
    for v in range(-maxFreqH, maxFreqH + 1):
        entryH = v + maxFreqH 
        for x in range(0, width):
            for y in range(0, height):
                 
                reconstruction[y,x] += (coeff[entryH, entryW][0] / (m*n)) * (cos(x * ww * u) * cos(y * wh * v) - sin(x * ww * u) * sin(y * wh * v)) - \
                                       (coeff[entryH, entryW][1] / (m*n)) * (cos(x * ww * u) * sin(y * wh * v) + sin(x * ww * u) * cos(y * wh * v))   


showImageF(reconstruction)
     
# Power
power = createImageF( 1 + 2 * maxFreqW, 1 + 2 * maxFreqH)
for kw,kh in itertools.product(range(-maxFreqW, maxFreqW + 1), range(-maxFreqH, maxFreqH + 1)):                  
    entryW = kw + maxFreqW 
    entryH = kh + maxFreqH 
    power[entryH, entryW] = sqrt(coeff[entryH, entryW][0] * coeff[entryH, entryW][0] +              \
                                 coeff[entryH, entryW][1] * coeff[entryH, entryW][1])
                 
    power[entryH, entryW] = log(1.0 + power[entryH, entryW])

# Show the log of the power 
powerLog = imageLogF(power)
showImageF(powerLog)

# Phase
phase = createImageF( 1 + 2 * maxFreqW, 1 + 2 * maxFreqH)
for kw,kh in itertools.product(range(-maxFreqW, maxFreqW + 1), range(-maxFreqH, maxFreqH + 1)):
    indexInArrayW = kw + maxFreqW 
    indexInArrayH = kh + maxFreqH 
    phase[indexInArrayH, indexInArrayW] = atan2(coeff[indexInArrayH, indexInArrayW][1],             \
                                                coeff[indexInArrayH, indexInArrayW][0])

# Plot phase
showImageF(phase)
