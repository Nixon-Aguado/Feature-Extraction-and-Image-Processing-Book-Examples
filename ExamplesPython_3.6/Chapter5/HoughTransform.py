'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 5
HoughTransform: Line detection by the Hough transform 
'''

# Set module functions
from ImageUtilities import imageReadL, createImageF, showImageF, showImageL, createScaleImageL
from ImageOperatorsUtilities import applyCannyEdgeDetector
from ImagePropertiesUtilities import imageMaxMin, peakDetectorImageL 
from PlotUtilities import plot3DHistogram 

# Math and iteration
from math import pi, tan
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    gaussianKernelSize = Gaussian kernel size. Filter noise
    sobelKernelSize = Sobel kernel size. Edge detection
    upperT = Upper threshold
    lowerT = Lower threshold
    peakDetection = Percentage of the maximum peak value that is considered for threshold
'''
pathToDir = "../../Images/Chapter5/Input/"
imageName = "Road.png"
gaussianKernelSize = 7
sobelKernelSize = 3
upperT = 0.5
lowerT = 0.3
peakDetection = 0.7

# Read image into array and show
inputImage, width, height = imageReadL(pathToDir + imageName)
showImageL(inputImage)

# Compute edges
magnitude, angle = applyCannyEdgeDetector(inputImage, gaussianKernelSize, sobelKernelSize, upperT, lowerT)
showImageF(magnitude)

# Two accumulators, for horizontal and vertical lines. Each one stores a range of 90 degrees
# The intersection c corresponds to the intersections with the lines x=0 and y=0 
accHorizontal = createImageF(2*height,90)
accVertical = createImageF(2*width,90);

# Gather evidence 
for x,y in itertools.product(range(0, width), range(0, height)):
    if magnitude[y,x] != 0:
        for m in range(0,90):
            
            # Lines between -45 and 45 degrees
            angle = ((-45 + m) * pi) / 180.0
            c = y - tan(angle) * x
            bucket = int(c)
            if bucket> 0 and bucket < 2*height - 1:
                weight = c - int(c)
                accHorizontal[m, bucket] += (1.0 - weight)
                accHorizontal[m, bucket+1] += weight
            
            # Lines between 45 and 135 degrees
            angle = ((45.0 + m) * pi) / 180.0
            c = x - y / tan(angle)
            bucket = int(c)
            if bucket> 0 and bucket < 2*width - 1:
                weight = c - int(c)
                accVertical[m, bucket] += (1.0 - weight)
                accVertical[m, bucket+1] += weight
            
# Find maximum
maxH, _ = imageMaxMin(accHorizontal)
maxV, _ = imageMaxMin(accVertical)
maximum = max(maxH, maxV)
peakThreshold = peakDetection * maximum

# Plot accumulators
plot3DHistogram(accHorizontal, [0,maximum])
plot3DHistogram(accVertical, [0,maximum])

# Prepare output image as a dark version of the input
outputImage = createScaleImageL(inputImage, 0.5)
    
# Peak detection
peakHorizontal = peakDetectorImageL(accHorizontal, peakThreshold)
peakVertical = peakDetectorImageL(accVertical, peakThreshold)

# Draw lines on output image
for peakIndex in range(0,len(peakHorizontal)):
    m = (peakHorizontal[peakIndex])[0]
    c = (peakHorizontal[peakIndex])[1]     
    strength = int(255.0 * accHorizontal[m, c] / maximum)
    angle = ((-45 + m) * pi) / 180.0
    for x in range(0, width -1):
        y = int(c + tan(angle) * x)
        if y > 0 and y < height -1:
            outputImage[y,x] = strength
            outputImage[y+1,x] = strength
    
for peakIndex in range(0,len(peakVertical)):
    m = (peakVertical[peakIndex])[0]
    c = (peakVertical[peakIndex])[1]  
        
    strength = int(255.0 * accVertical[m, c] / maximum)
    angle = ((45 + m) * pi) / 180.0
    for y in range(0, height -1):
        x = int(c + y / tan(angle))
        if x > 0 and x < width -1:
            outputImage[y,x] = strength
            outputImage[y,x+1] = strength  
                
showImageL(outputImage)

