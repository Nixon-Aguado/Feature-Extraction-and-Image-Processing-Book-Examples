'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
CurveatureDetection: Obtain curvature by computing angle differences or
'''

# Set module functions
from ImageUtilities import imageReadL, createImageF, showImageF, showImageL
from ImageOperatorsUtilities import applyCannyEdgeDetector

# Math and iteration
from math import cos, sin
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    GaussianKernelSize = Gaussian kernel size. Filter noise
    sobelKernelSize = Sobel kernel size. Edge detection
    upperT = Upper threshold
    lowerT = Lower threshold
    windowDelta = Size of the region used to find neighbors
'''
pathToDir = "../../Images/Chapter4/Input/"
imageName = "Shapes.png"
GaussianKernelSize = 7
sobelKernelSize = 3
upperT = 0.4
lowerT = 0.2
windowDelta = 2

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Compute edges
magnitude, angle = applyCannyEdgeDetector(inputImage, GaussianKernelSize, sobelKernelSize, upperT, lowerT)
showImageF(magnitude)

# Compute curvature by subtracting the direction of neighbors
curvature = createImageF(width, height)
for x,y in itertools.product(range(0, width), range(0, height)):
    # Edge
    if magnitude[y,x] > 0:  
        # Consider neighbor edges
        edgesNeigbor = [ ]
        for wx,wy in itertools.product(range(-windowDelta, windowDelta+1),       \
                                       range(-windowDelta, windowDelta+1)):
            if magnitude[y+wy, x+wx] > 0 : 
                edgesNeigbor.append((y+wy,x+wx))
               
        # Use dot product to measure angle difference
        np = len(edgesNeigbor)     
        for p in range(0, np):
            y1 = (edgesNeigbor[p])[0]
            x1 = (edgesNeigbor[p])[1]
            curvature[y,x] += 1.0-(cos(angle[y1,x1]) * cos(angle[y,x])         \
                                 + sin(angle[y1,x1]) * sin(angle[y,x]))
        if np > 0:
            curvature[y,x] /= np 

showImageF(curvature, 2)
