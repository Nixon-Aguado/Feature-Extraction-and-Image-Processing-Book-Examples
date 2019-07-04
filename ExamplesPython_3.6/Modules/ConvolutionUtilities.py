'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
ConvolutionUtilities: Helper module to apply convolutions. Functions to create kernels and perform convolutions
'''

# Images
from ImageUtilities import createImageF, createImageL
from ImagePropertiesUtilities import imageMaxMin

# Math and iteration
from math import exp, pow, factorial, sqrt, atan2
from timeit import itertools

# Generate a Gaussian kernel
def createGaussianKernel(kernelSize):
    
    sigma = kernelSize / 3.0
    
    kernelImage = createImageF(kernelSize, kernelSize)
    centre = (kernelSize - 1) / 2
    for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
        kernelImage[y,x] =  exp( -0.5 * (pow((x - centre)/sigma, 2.0) + \
                                         pow((y - centre)/sigma, 2.0)) ) 
        
    return kernelImage

# Create a Sobel kernel k
def createSobelKernel(kenrelSize):
    
    sobelX = createImageF(kenrelSize, kenrelSize)
    sobelY = createImageF(kenrelSize, kenrelSize)
    
    # Create kernel
    for x,y in itertools.product(range(0, kenrelSize), range(0, kenrelSize)):
        
        # Smooth
        smoothX = factorial(kenrelSize - 1) / (factorial(kenrelSize - 1 - x) * factorial(x))
        smoothY = factorial(kenrelSize - 1) / (factorial(kenrelSize - 1 - y) * factorial(y))
    
        # Pascal 
        if ( kenrelSize - 2 - x >= 0):
            p1X = factorial(kenrelSize - 2) / (factorial(kenrelSize - 2 - x) * factorial(x))
        else:
            p1X = 0
        
        if ( kenrelSize - 2 - y >= 0):    
            p1Y = factorial(kenrelSize - 2) / (factorial(kenrelSize - 2 - y) * factorial(y))
        else:
            p1Y = 0
            
        # Pascal shift to the right
        xp = x - 1
        if ( kenrelSize - 2 - xp >= 0 and xp >= 0):
            p2X = factorial(kenrelSize - 2) / (factorial(kenrelSize - 2 - xp) * factorial(xp))
        else:
            p2X = 0
            
        yp = y - 1
        if ( kenrelSize - 2 - yp >= 0 and yp >= 0):
            p2Y = factorial(kenrelSize - 2) / (factorial(kenrelSize - 2 - yp) * factorial(yp))
        else:
            p2Y = 0
            
        # Sobel    
        sobelX[y,x] = smoothX * (p1Y - p2Y)
        sobelY[y,x] = smoothY * (p1X - p2X)                                        

    return sobelX, sobelY

# Create a Laplacian kernel 
def createLaplacianKernel(kernelSize, sigma):

    # kernel
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
                                                                               
    return kernelLaplacian

# Apply kernel to an image returning a gray level image
def applyKernel(inputImage, kernelImage):
    height = len(inputImage)
    width = len(inputImage[0])    

    kernelHeight = len(kernelImage)
    kerelWidth = len(kernelImage[0]) 
    
    kernelCentreY = int((kernelHeight - 1) / 2)
    kernelCentreX = int((kerelWidth - 1) / 2)
    
    # Create images to store the result
    outputImage = createImageL(width, height)
    
    for x,y in itertools.product(range(0, width), range(0, height)):
        sumKernel = 0
        sumKernelWeights = 0
        for wx,wy in itertools.product(range(0, kerelWidth), range(0, kernelHeight)):
            posY = y + wy - kernelCentreY
            posX = x + wx - kernelCentreX 
            
            if posY > -1 and posY <  height and  posX > -1 and posX <  width:
                sumKernel += float(inputImage[posY,posX]) * kernelImage[wy, wx]
                sumKernelWeights += kernelImage[wy, wx]
        
        if sumKernelWeights > 0:
            outputImage[y,x] = sumKernel / sumKernelWeights
                
    return outputImage

# Apply kernel returning a float image
def applyKernelF(inputImage, kernelImage):
    height = len(inputImage)
    width = len(inputImage[0])    

    kernelHeight = len(kernelImage)
    kernelWidth = len(kernelImage[0]) 
    
    kernelCentreY = int((kernelHeight - 1) / 2)
    kernelCentreX = int((kernelWidth - 1) / 2)
    
    # Create images to store the result
    outputImage = createImageF(width, height)
    
    for x,y in itertools.product(range(0, width), range(0, height)):
        sumKernel = 0.0
        sumKernelWeights = 0.0
        for wx,wy in itertools.product(range(0, kernelWidth), range(0, kernelHeight)):
            posY = y + wy - kernelCentreY
            posX = x + wx - kernelCentreX 
            
            if posY > -1 and posY <  height and  posX > -1 and posX <  width:
                sumKernel += float(inputImage[posY,posX]) * float(kernelImage[wy, wx])
                sumKernelWeights += float(kernelImage[wy, wx])
        
        # If we have to normalize
        if sumKernelWeights != 0.0:
            outputImage[y,x] = sumKernel / sumKernelWeights
        else:
            outputImage[y,x] = sumKernel
            
    return outputImage

# Apply kernels to an image return magnitude and angle
def applyKernelMA(inputImage, kernelX, kernelY, normalizeMagnitude = False):
    height = len(inputImage)
    width = len(inputImage[0])    

    kernelHeight = len(kernelX)
    kerelWidth = len(kernelX[0]) 
    
    kernelCentreY = int((kernelHeight - 1) / 2)
    kernelCentreX = int((kerelWidth - 1) / 2)
    
    # Create images to store the result
    magnitude = createImageF(width, height)
    direction = createImageF(width, height)
    imageX = createImageF(width, height)
    imageY = createImageF(width, height)
 
    # Convolution with two kernels
    for x,y in itertools.product(range(kernelCentreX, width - kernelCentreX),          \
                                 range(kernelCentreY, height - kernelCentreY)):    
        sumKernel = [0.0, 0.0]
        sumKernelWeights = [0.0, 0.0]
        for wx,wy in itertools.product(range(0, kerelWidth), range(0, kernelHeight)):

            posY = y + wy - kernelCentreY
            posX = x + wx - kernelCentreX 
                
            if posY > -1 and posY <  height and  posX > -1 and posX <  width:
                sumKernel[0] += float(inputImage[posY,posX]) * float(kernelX[wy, wx])
                sumKernelWeights[0] += float(kernelX[wy, wx])
                
                sumKernel[1] += float(inputImage[posY,posX]) * float(kernelY[wy, wx])
                sumKernelWeights[1] += float(kernelY[wy, wx])
            
        # If we have to normalize
        if sumKernelWeights[0] != 0.0:
            imageX[y,x] = sumKernel[0] / sumKernelWeights[0]
        else:
            imageX[y,x] = sumKernel[0]
            
        # If we have to normalize
        if sumKernelWeights[1] != 0.0:
            imageY[y,x] = sumKernel[1] / sumKernelWeights[1]
        else:
            imageY[y,x] = sumKernel[1]
                
        magnitude[y,x] = sqrt(imageX[y,x] * imageX[y,x] + imageY[y,x] * imageY[y,x])   
        direction[y,x] = atan2(imageY[y,x], imageX[y,x])   
            
    if normalizeMagnitude == True:
        maximum, minimum = imageMaxMin(magnitude)
        
        for x,y in itertools.product(range(0, width), range(0, height)):
            magnitude[y,x] = (magnitude[y,x] - minimum) / float(maximum - minimum)
    
    return magnitude, direction, imageX, imageY
