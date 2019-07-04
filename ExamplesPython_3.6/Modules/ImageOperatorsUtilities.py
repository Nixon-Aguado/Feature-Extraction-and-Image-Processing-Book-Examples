'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
ImageOperatorsUtilities: Helper module to process an image
'''

# Images
from ImageUtilities import createImageF, createImageL, createVectorI
from ConvolutionUtilities import createGaussianKernel, createSobelKernel, applyKernelF, applyKernelMA

# Array to store image data
from numpy import amax

# Math
from math import log, pi, sin, cos

# Iteration
from timeit import itertools

# Create a log of an array containing float data
def imageLogF(image, scale = 100):
    height = len(image)
    width = len(image[0])
    
    maximum = amax(image) 
    
    # Create a gray level image for display
    logImage = createImageF(width ,height)
    
    # Scale the float values into gray scale values without 0 and negative values
    for y in range(0, height):
        for x in range(0, width):
            # Set to log value 
            logImage[y, x] = log(1.0 + (abs(image[y, x]) / maximum) * scale)
  
    return logImage

# Compute the histogram of an image
def computeHistogram(inputImage):
    height = len(inputImage)
    width = len(inputImage[0])   
    
    # Vector of integers values to store the number of times a pixel value is repeated
    outputHistogram = createVectorI(256)
  
    # Get the number of times a pixel value is found in the image
    for x,y in itertools.product(range(0, width), range(0, height)):
        pixelValue = inputImage[y,x]
        outputHistogram[pixelValue] += 1
            
    return outputHistogram 

# Return byte thresholded image
def thresholdImage(inputImage, threshold, binary = True):

    height = len(inputImage)
    width = len(inputImage[0])   
    
    # Create images to store the result
    outputImage = createImageL(width, height)

    # Set the pixels in the output image according to the accumulate histogram
    for x,y in itertools.product(range(0, width), range(0, height)):
  
        if inputImage[y,x] > threshold:
            if binary:
                outputImage[y,x] = 255
            else:
                outputImage[y,x] = inputImage[y,x]
        else:
            outputImage[y,x] = 0
                
    return outputImage


# Apply Canny operator to an image
def applyCannyEdgeDetector(inputImage, GaussianKernelSize, sobelKernelSize, upperT, lowerT, returnGradient = False):
    
    height = len(inputImage)
    width = len(inputImage[0])    
    
    normalizeMagnitude = True
    windowDelta = 1
    
    # Apply Gaussian kernel
    gaussianKernel = createGaussianKernel(GaussianKernelSize)
    gaussianImage = applyKernelF(inputImage, gaussianKernel)
    
    # Apply Sobel kernel. We use normalized magnitude in this example
    sobelX, sobelY = createSobelKernel(sobelKernelSize)
    magnitude, angle, mX, mY = applyKernelMA(gaussianImage, sobelX, sobelY, normalizeMagnitude)
    
    # Weight magnitude by the variance. This is useful for corner extraction since suppress the internal corner
    weightedMagnitude = createImageF(width, height)
    for x,y in itertools.product(range(0, width), range(0, height)):
        sumKernel = 1.0/8.0
        for wx,wy in itertools.product(range(-1,2), range(-1, 2)):
 
            posY = y + wy
            posX = x + wx 
            
            if posY > -1 and posY <  height and  posX > -1 and posX <  width:
                sumKernel += abs(float(inputImage[posY,posX]) -  float(inputImage[y,x]))
    
        sumKernel /= 8.0
        weightedMagnitude[y,x] = magnitude[y,x] *  sumKernel

    # To store maximum suppression image
    maxImage = createImageF(width, height)
    
    # Non-maximum suppression
    border = GaussianKernelSize
    for x,y in itertools.product(range(border, width - border), 
                                 range(border, height - border)):
        
        # Only potential edges can be maximum
        if magnitude[y,x] > lowerT:     
           
            # The normal angle is perpendicular to the edge angle
            normalAngle = angle[y,x] - pi / 2.0
            
            # Make sure the angle is between 0 and pi 
            while normalAngle < 0:
                normalAngle += pi
            while normalAngle > pi:
                normalAngle -= pi
            
            # Angle defining the first point
            baseAngle = int( 4 * normalAngle / pi ) * (pi / 4.0)
            
            # Integer delta positions for interpolation
            # We use -y since the image origin is in top corner
            x1, y1 = int(round(cos(baseAngle))), -int(round(sin(baseAngle)))
            x2, y2 = int(round(cos(baseAngle + pi / 4.0))),                     \
                    -int(round(sin(baseAngle + pi / 4.0)))
            
            # How far we are from (x1,y1). Maximum difference is math.pi / 4.0, so we multiply by 2
            w = cos(2.0*(normalAngle - baseAngle))
            
            # Point to interpolate
            M1 = w * weightedMagnitude[y+y1,x+x1] + (1.0 - w) * weightedMagnitude[y+y2,x+x2]
            
            # Point to interpolate for pixels in the other side of the edge
            M2 = w * weightedMagnitude[y-y1,x-x1] + (1.0 - w) * weightedMagnitude[y-y2,x-x2]
             
            # Determine if it is a maximum. If so make sure it will be preserved 
            if weightedMagnitude[y,x] > M1 and weightedMagnitude[y,x] > M2:
                maxImage[y,x] = magnitude[y,x]

    # To compute hysteresis thresholded images we require two thresholds
    edges = createImageF(width, height)
    potentialEdges = [ ]

    # Divide pixels as edges, no edges and we are not sure
    for x,y in itertools.product(range(1, width-1), range(1, height-1)):            
        # These are edges
        if maxImage[y,x] > upperT:
            edges[y,x] = 255
            
        # These are pixels that we do not want as edges   
        if maxImage[y,x] < lowerT:
            edges[y,x] = 0
            
        # These may be edges    
        if maxImage[y,x] > lowerT and maxImage[y,x] <= upperT:
            edges[y,x] = 128

    # Resolve the potential edges
    for x,y in itertools.product(range(1, width-1), range(1, height-1)):            
        # For each edge
        if edges[y,x] == 255:
            
            # Examine neighbour
            potentialEdges = [ ]
            for wx,wy in itertools.product(range(-windowDelta, windowDelta+1), range(-windowDelta, windowDelta+1)):                  
                # It becomes an edge
                if edges[y+wy,x+wx] == 128:
                    edges[y+wy,x+wx] = 255
                    potentialEdges.append((y+wy,x+wx))
            
            # Look into new edges
            while len(potentialEdges) > 0:
                # Take element from potential edges
                y = (potentialEdges[0])[0]
                x = (potentialEdges[0])[1]
                potentialEdges = potentialEdges[1:]
                
                # Examine neighbour
                for wx,wy in itertools.product(range(-windowDelta, windowDelta+1), range(-windowDelta, windowDelta+1)): 
                    # It becomes an edge
                    if edges[y+wy,x+wx] == 128:
                        edges[y+wy,x+wx] = 255
                        potentialEdges.append((y+wy,x+wx))
        
    # Clean up remaining potential edges                            
    for x,y in itertools.product(range(1, width-1), range(1, height-1)): 
        if edges[y,x] == 128:
            edges[y,x] = 0
    
    if returnGradient == False:
        return edges, angle 
 
    return edges, angle , mX, mY

