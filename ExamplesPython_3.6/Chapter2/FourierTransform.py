'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 2
FourierTransform: Compute the Fourier transform of an image and display the magnitude and phase
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from ImageOperatorsUtilities import imageLogF
from PrintUtilities import printProgress

# Math and iteration functions
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
maxFrequencyW = int(width /2)
maxFrequencyH = int(height/2)
numCoefficientsW = 1 + 2 * maxFrequencyW
numCoefficientsH = 1 + 2 * maxFrequencyH
coefficients = createImageF(numCoefficientsW ,numCoefficientsH , 2)

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

# Compute coefficients
for kw in range(-maxFrequencyW, maxFrequencyW + 1):
    printProgress(kw + maxFrequencyW, numCoefficientsW)
    indexInArrayW = kw + maxFrequencyW 
    for kh in range(-maxFrequencyH, maxFrequencyH + 1):
        indexInArrayH = kh + maxFrequencyH 
        for x,y in itertools.product(range(0, width), range(0, height)):
            coefficients[indexInArrayH, indexInArrayW][0] +=  inputImage[y,x] *          \
                 (cos(x * ww * kw) * cos(y * wh * kh) -  sin(x * ww * kw) * sin(y * wh * kh))                 
            coefficients[indexInArrayH, indexInArrayW][1] +=  inputImage[y,x] *           \
                 (cos(x * ww * kw) * sin(y * wh * kh) +  sin(x * ww * kw) * cos(y * wh * kh)) 
            
for kw in range(-maxFrequencyW, maxFrequencyW + 1):
    printProgress(kw + maxFrequencyW, numCoefficientsW)
    indexInArrayW = kw + maxFrequencyW 
    for kh in range(-maxFrequencyH, maxFrequencyH + 1):
        indexInArrayH = kh + maxFrequencyH 
        coefficients[indexInArrayH, indexInArrayW][0] *= m*n 
        coefficients[indexInArrayH, indexInArrayW][1] *= m*n
                            
# Power
power = createImageF( 1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH)
for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),                  \
                               range(-maxFrequencyH, maxFrequencyH + 1)):                  
    indexInArrayW = kw + maxFrequencyW 
    indexInArrayH = kh + maxFrequencyH 
    power[indexInArrayH, indexInArrayW] =                                                 \
            sqrt(coefficients[indexInArrayH, indexInArrayW][0] *                          \
                 coefficients[indexInArrayH, indexInArrayW][0] +                          \
                 coefficients[indexInArrayH, indexInArrayW][1] *                          \
                 coefficients[indexInArrayH, indexInArrayW][1])
            
# Show the log of the power 
powerLog = imageLogF(power)
showImageF(powerLog)

# Phase
phase = createImageF( 1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH)
for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),                  \
                               range(-maxFrequencyH, maxFrequencyH + 1)):
    indexInArrayW = kw + maxFrequencyW 
    indexInArrayH = kh + maxFrequencyH 
    phase[indexInArrayH, indexInArrayW] =                                                 \
                    atan2(coefficients[indexInArrayH, indexInArrayW][1],                  \
                          coefficients[indexInArrayH, indexInArrayW][0])

# Show phase
showImageF(phase)
