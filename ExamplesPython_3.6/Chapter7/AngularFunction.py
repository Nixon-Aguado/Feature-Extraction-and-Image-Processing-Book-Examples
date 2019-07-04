'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 7
AngularFunction: Compute the angular functions of a shape
'''
        
# Set module functions
from ImageUtilities import imageReadL, createImageF, showImageF, showImageL
from ImageRegionsUtilities import findLongestCentredSegmentinImage, showShapeinImage
from PlotUtilities import plotCurveXY 

# Math and iteration
from math import pi, sqrt, atan2
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    gaussianKernelSize = Gaussian kernel size. Filter noise
    sobelKernelSize = Sobel kernel size. Edge detection
    upperT = Upper threshold
    lowerT = Lower threshold
'''
pathToDir = "../../Images/Chapter7/Input/"
imageName = "Shape.png"
gaussianKernelSize = 5
sobelKernelSize = 3
upperT = 0.3
lowerT = 0.05

# Obtain a shape from the input image and draw it
centre, shape, width, height = findLongestCentredSegmentinImage(pathToDir + imageName,   \
                                        gaussianKernelSize, sobelKernelSize, upperT, lowerT)
showShapeinImage(shape, centre, width, height)

# Compute the accumulative arc lengths 
numPoints = len(shape[0])
sumLenghts = []  
y0, x0 = shape[0, numPoints-1], shape[1, numPoints-1]
shapeLenght = 0.0
for p in range(0, numPoints):
    y,x = shape[0,p], shape[1,p]
    shapeLenght += sqrt((y-y0)*(y-y0) + (x-x0)*(x-x0))
    sumLenghts.append(shapeLenght)
    y0,x0 = y,x
  
# Normalised arc lengths
normLenghts = []
for p in range(0, numPoints):
    normLenghts.append((2.0*pi*sumLenghts[p])/shapeLenght);

# Compute angular function by an average window
windowSize = [1,10]
d = float(windowSize[1] -windowSize[0])  
angularFunc = [ ]
for p in range(0, numPoints):
    x1,x2,y1,y2 = 0.0, 0.0, 0.0, 0.0
    # Average change 
    for q in range(windowSize[0], windowSize[1]):
        pa,pb = p-q,p+q
        if pa<0:           pa += numPoints
        if pb>=numPoints:  pb -= numPoints
        
        ya,xa = shape[0,pa], shape[1,pa]
        yb,xb = shape[0,pb], shape[1,pb]
        
        x1,y1 = x1+xa, y1+ya
        x2,y2 = x2+xb, y2+yb    
    dx, dy = (x2-x1)/d, (y2-y1)/d
    angle = atan2(dy, dx)
    angularFunc.append(angle)

# Compute cumulative angular function
cumulativeFunc = [ ]
angle0 = angularFunc[numPoints-1]
sumAngle = 0.0
for p in range(0, numPoints):
    angle = angularFunc[p]
    if abs(angle-angle0) < pi:
        sumAngle += angle-angle0
    else:
        sumAngle += angle-(angle0 + 2.0 *pi)   
    cumulativeFunc.append(sumAngle)
    angle0 = angle    

# Compute cumulative angular accumulated
cumNormFunc = [ ]
for p in range(0, numPoints):
    cumNormFunc.append(cumulativeFunc[p]+normLenghts[p])
    
plotCurveXY(sumLenghts,angularFunc, [-3.2 , 3.2])
plotCurveXY(sumLenghts,cumulativeFunc, [-7, 0])
plotCurveXY(sumLenghts,cumNormFunc, [-3, 3])
