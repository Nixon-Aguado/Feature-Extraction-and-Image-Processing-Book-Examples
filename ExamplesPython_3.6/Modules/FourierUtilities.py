'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
FourierUtilties: Helper module for Fourier analysis
'''

# Images
from ImageUtilities import createImageF
from PrintUtilities import printProgress

# Math and iteration functions
from math import pi, sin, cos, atan2, sqrt
from timeit import itertools

# Compute the power and phase from an input image
def computePowerandPhase(inputImage):
    height = len(inputImage)
    width = len(inputImage[0])   
    
    # Create coefficients Image. Two floats to represent a complex number
    # Maximum frequency according to sampling
    maxFrequencyW = int(width /2)
    maxFrequencyH = int(height/2)
    numCoefficientsW = 1 + 2 * maxFrequencyW
    numCoefficientsH = 1 + 2 * maxFrequencyH
    coefficients = createImageF(numCoefficientsW ,numCoefficientsH , 2)
    
    # Adjust the size of the data to be even
    m = float(width)
    n = float(height)
    if width % 2 == 0:
        m = width + 1
    if height % 2 == 0:
        n = height + 1  
    
    # Fundamental frequency
    ww = (2.0 * pi) / float(m)
    wh = (2.0 * pi) / float(n)
    
    # Compute values
    for kw in range(-maxFrequencyW, maxFrequencyW + 1):
        printProgress(kw + maxFrequencyW, numCoefficientsW)
        for kh in range(-maxFrequencyH, maxFrequencyH + 1):
            indexInArrayW = kw + maxFrequencyW 
            indexInArrayH = kh + maxFrequencyH 

            for x,y in itertools.product(range(0, width), range(0, height)):
                coefficients[indexInArrayH, indexInArrayW][0] += inputImage[y,x] * (cos(x * ww * kw) * cos(y * wh * kh) - sin(x * ww * kw) * sin(y * wh * kh))
                coefficients[indexInArrayH, indexInArrayW][1] += inputImage[y,x] * (cos(x * ww * kw) * sin(y * wh * kh) + sin(x * ww * kw) * cos(y * wh * kh))
    
    # Power
    powerImage = createImageF( 1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH) 
    
    for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),                  \
                               range(-maxFrequencyH, maxFrequencyH + 1)):   
        indexInArrayW = kw + maxFrequencyW 
        indexInArrayH = kh + maxFrequencyH 
        powerImage[indexInArrayH, indexInArrayW] = sqrt(coefficients[indexInArrayH, indexInArrayW][0] * coefficients[indexInArrayH, indexInArrayW][0] + \
                                                        coefficients[indexInArrayH, indexInArrayW][1] * coefficients[indexInArrayH, indexInArrayW][1])
    
    # Phase
    phaseImage = createImageF( 1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH)
    
    for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),             \
                                   range(-maxFrequencyH, maxFrequencyH + 1)):    
        indexInArrayW = kw + maxFrequencyW 
        indexInArrayH = kh + maxFrequencyH 
        phaseImage[indexInArrayH, indexInArrayW] = atan2(coefficients[indexInArrayH, indexInArrayW][1], coefficients[indexInArrayH, indexInArrayW][0])
    
    return powerImage, phaseImage 

# Compute Fourier transform coefficients from an image
def computeCoefficients(inputImage):
    height = len(inputImage)
    width = len(inputImage[0])   
    
    # Create coefficients Image. Two floats to represent a complex number
    # Maximum frequency according to sampling
    maxFrequencyW = int(width /2)
    maxFrequencyH = int(height/2)
    numCoefficientsW = 1 + 2 * maxFrequencyW
    numCoefficientsH = 1 + 2 * maxFrequencyH
    coefficients = createImageF(numCoefficientsW ,numCoefficientsH , 2)

    # Adjust the size of the data to be even
    m = float(width)
    n = float(height)
    if width % 2 == 0:
        m = width + 1
    if height % 2 == 0:
        n = height + 1
    
    # Fundamental frequency
    ww = (2.0 * pi) / float(m)
    wh = (2.0 * pi) / float(n)
    
    # Compute values
    for kw in range(-maxFrequencyW, maxFrequencyW + 1):
        printProgress(kw + maxFrequencyW, numCoefficientsW - 1)
        for kh in range(-maxFrequencyH, maxFrequencyH + 1):
            indexInArrayW = kw + maxFrequencyW 
            indexInArrayH = kh + maxFrequencyH 
            
            for x,y in itertools.product(range(0, width), range(0, height)):
                coefficients[indexInArrayH, indexInArrayW][0] += inputImage[y,x] * (cos(x * ww * kw) * cos(y * wh * kh) - sin(x * ww * kw) * sin(y * wh * kh))
                coefficients[indexInArrayH, indexInArrayW][1] += inputImage[y,x] * (cos(x * ww * kw) * sin(y * wh * kh) + sin(x * ww * kw) * cos(y * wh * kh))
  
    for kw in range(-maxFrequencyW, maxFrequencyW + 1):
        for kh in range(-maxFrequencyH, maxFrequencyH + 1):
            coefficients[indexInArrayH, indexInArrayW][0] /= (m*n)
            coefficients[indexInArrayH, indexInArrayW][1] /= (m*n)
  
    return coefficients, maxFrequencyW, maxFrequencyH 

# Return the power image from the coefficients
def computePowerfromCoefficients(coefficients):
    # Maximum frequency
    maxFrequencyH = int((len(coefficients) - 1) / 2)
    maxFrequencyW = int((len(coefficients[0]) - 1) / 2) 
      
    # Power
    powerImage = createImageF( 1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH) 
    
    for kw in range(-maxFrequencyW, maxFrequencyW + 1):
        printProgress(kw + maxFrequencyW, 2 * maxFrequencyW)
        for kh in range(-maxFrequencyH, maxFrequencyH + 1):
            indexInArrayW = kw + maxFrequencyW 
            indexInArrayH = kh + maxFrequencyH 
            powerImage[indexInArrayH, indexInArrayW] = sqrt(coefficients[indexInArrayH, indexInArrayW][0] * coefficients[indexInArrayH, indexInArrayW][0] + \
                                                            coefficients[indexInArrayH, indexInArrayW][1] * coefficients[indexInArrayH, indexInArrayW][1])
     
    return powerImage 

# Inverse transform
def reconstruction(coefficients): 
    # Maximum frequency
    maxFrequencyH = int((len(coefficients) - 1) / 2)
    maxFrequencyW = int((len(coefficients[0]) - 1) / 2) 
    
    height = 2 * maxFrequencyH
    width = 2 * maxFrequencyW 
        
    # Adjust the size of the data to be even
    m = float(width)
    n = float(height)
    if width % 2 == 0:
        m = width + 1
    if height % 2 == 0:
        n = height + 1  
        
    # Fundamental frequency
    ww = (2.0 * pi) / float(m)
    wh = (2.0 * pi) / float(n)

    reconstructionImage = createImageF(m, n)
    for x in range(0, width):
        printProgress(x, width - 1)
        for y in range(0, height):
            for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),             \
                                           range(-maxFrequencyH, maxFrequencyH + 1)):    
                indexInArrayW = kw + maxFrequencyW 
                indexInArrayH = kh + maxFrequencyH 
                reconstructionImage[y,x] += \
                        (coefficients[indexInArrayH, indexInArrayW][0] / (m*n)) * (cos(x * ww * kw) * cos(y * wh * kh) - sin(x * ww * kw) * sin(y * wh * kh)) + \
                        (coefficients[indexInArrayH, indexInArrayW][1] / (m*n)) * (cos(x * ww * kw) * sin(y * wh * kh) + sin(x * ww * kw) * cos(y * wh * kh))   

    return reconstructionImage

