'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 2
LowHighFilters: Filter an image in frequency domain 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from FourierUtilities import computeCoefficients, computePowerfromCoefficients, reconstruction
from ImageOperatorsUtilities import imageLogF

# Math functions and iteration
from math import sqrt
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter2/Input/"
imageName = "Giraffe.png" 

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show images
showImageL(inputImage)

# Compute coefficients 
coeff, maxFreqW, maxFreqH = computeCoefficients(inputImage)

# Create low and high versions
coeffLow = createImageF( 1 + 2 * maxFreqW, 1 + 2 * maxFreqH, 2)
coeffHigh = createImageF( 1 + 2 * maxFreqW, 1 + 2 * maxFreqH, 2)

# Filter
cutFrequency = maxFreqW / 8

for kw,kh in itertools.product(range(-maxFreqW, maxFreqW + 1),        \
                               range(-maxFreqH, maxFreqH + 1)):   
    IndexW, indexH = kw + maxFreqW, kh + maxFreqH
    
    if sqrt(kw * kw + kh * kh) < cutFrequency:
        coeffLow[indexH, IndexW][0] = coeff[indexH, IndexW][0]
        coeffLow[indexH, IndexW][1] = coeff[indexH, IndexW][1]
    else:
        coeffHigh[indexH, IndexW][0] = coeff[indexH, IndexW][0]
        coeffHigh[indexH, IndexW][1] = coeff[indexH, IndexW][1]
            
# Power
powerLow = computePowerfromCoefficients(coeffLow)
powerHigh = computePowerfromCoefficients(coeffHigh)

# Show power
powerLowLog = imageLogF(powerLow)
powerHighLog = imageLogF(powerHigh)
showImageF(powerLowLog)
showImageF(powerHighLog)

# Reconstruct image
imageLow = reconstruction(coeffLow)
imageHigh = reconstruction(coeffHigh)

showImageF(imageLow)
showImageF(imageHigh)

