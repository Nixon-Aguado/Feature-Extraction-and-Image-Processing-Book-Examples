'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 8
WaterShedEdgeTransform: Compute Watershed transform by considering the edge image
'''
        
# Set module functions
from ImageUtilities import imageReadL, showImageL,createImageF, showImageF
from ImageOperatorsUtilities import applyCannyEdgeDetector
from ImageRegionsUtilities import watherShed
from PrintUtilities import printProgress
 
# Math and iteration
from math import  sqrt
from _testcapi import FLT_MAX
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    kernelSize = Gaussian and Sobel kernel size
    normalizeMagnitude = Normalise the convolution output
    upperT = upper threshold
    lowerT = lower threshold
    windowDelta = Size of window used in hysteresis
    suppWindow = Size of the window used to find maxima    
'''
pathToDir = "../../Images/Chapter8/Input/"
imageName = "Logs.png"
cannyKernelSize = 7
upperT = 0.5
lowerT = 0.1
windowDelta = 3
suppWindow = 5

# Read image into array and show
inputImage, width, height = imageReadL(pathToDir+imageName)
showImageL(inputImage)

# Compute edges
magnitude, angle = applyCannyEdgeDetector(inputImage, cannyKernelSize, cannyKernelSize, upperT, lowerT)
showImageF(magnitude)

# Divide pixels into edge and region pixels
edgePixels = [ ]
shapeImage = [ ]
for x,y in itertools.product(range(0, width), range(0, height)):
    if magnitude[y,x] > 0:
        edgePixels.append((y,x))
    shapeImage.append((y,x))

# Radial is the minimal distance to the edge
distanceImage = createImageF(width, height)
numEdges = len(edgePixels)
for x in range(0, width):
    printProgress(x, width)
    for y in range(0, height):
        minEdgeDist = FLT_MAX
        for indexEdge in range(0, numEdges):
            edgeY, edgeX = (edgePixels[indexEdge])[0], (edgePixels[indexEdge])[1]
            minEdgeDist = min(minEdgeDist, sqrt((edgeX-x)**2+(edgeY-y)**2) )
        # We define an edge in a distance image as 1. In this case we do not have edges so all flood
        distanceImage[y,x] = minEdgeDist + 2.0
showImageF(distanceImage)

# Watershed of the distance image
watershed = watherShed(distanceImage, shapeImage, suppWindow)
showImageF(watershed)     

