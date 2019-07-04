'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 5
InvariantGeneralizedHoughTransform: Shape detection by the invariant generalized Hough transform 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageF, showImageL, createImageF
from ImageOperatorsUtilities import applyCannyEdgeDetector
from ImageRegionsUtilities import computeReferencePoint
from PlotUtilities import plot3DHistogram 

# Math and iteration
from math import pi, tan, sqrt, atan
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    gaussianKernelSize = Gaussian kernel size. Filter noise
    sobelKernelSize = Sobel kernel size. Edge detection
    upperT = Upper threshold
    lowerT = Lower threshold
    numEntries = Size of the R table
    minimaDistPoints = To avoid to use close pairs of points
    maxDistPoints = To gather evidence in local regions
'''
pathToDir = "../../Images/Chapter5/Input/"
imageName = "Works.png"
templateName = "TemplateWorks.png"
gaussianKernelSize = 5
sobelKernelSize = 3
upperT = 0.4
lowerT = 0.3
numEntries = 90
minimaDistPoints = 20
maxDistPoints = 200

# Defines the direction of the line used to find pairs of points
alpha = pi/2.0

# Read image into array and show
templateImage, widthTemplate, heightTemplate = imageReadL(pathToDir + templateName)
inputImage, width, height = imageReadL(pathToDir + imageName)
showImageL(templateImage)
showImageL(inputImage)

# Compute edges
magnitudeTemplate, angleTemplate = applyCannyEdgeDetector(templateImage, gaussianKernelSize, sobelKernelSize, upperT, lowerT)
magnitude, angle = applyCannyEdgeDetector(inputImage, gaussianKernelSize, sobelKernelSize, upperT, lowerT)
showImageF(magnitudeTemplate)
showImageF(magnitude)

# Compute reference point in the template. Template centre
refPoint, edgePoints = computeReferencePoint(magnitudeTemplate)

# Find the pairs of points in the template according to alpha angle
pairPoints = []
numPts = len(edgePoints)
for p in range(0, numPts):
    y1, x1 = (edgePoints[p])[0], (edgePoints[p])[1]
    # We are looking for two points along the line with slope m
    m = tan(angleTemplate[y1,x1] - alpha)
    if m>-1 and m<1:
        if x1 < refPoint[1]:
            xi, xf, step = x1 + minimaDistPoints,widthTemplate, 1
        else:
            xi, xf, step = x1 - minimaDistPoints, 0, -1
        for x2 in range(xi, xf, step):
            y2 = int(y1 - m * (x2 - x1))
            if y2 > 0 and y2 < heightTemplate and magnitudeTemplate[y2,x2] != 0:
                pairPoints.append((y2,x2,y1,x1))
                break        
    else:
        m = 1.0/m
        if y1 < refPoint[0]:
            yi, yf, step = y1 + minimaDistPoints,heightTemplate, 1
        else:
            yi, yf, step = y1 - minimaDistPoints, 0, -1
        for y2 in range(yi, yf, step):
            x2 = int(x1 - m * (y2 - y1))
            if x2 > 0 and x2 < widthTemplate and magnitudeTemplate[y2,x2] != 0:
                pairPoints.append((y2,x2,y1,x1))
                break    
                 
# Build table (k,c) from each pair of points 
rTable = [[] for entryIndex in range(numEntries)]     
deltaAngle = pi / (numEntries - 1.0)  
numPairs = len(pairPoints)
for pair in range(0, numPairs):
    y2, x2 = (pairPoints[pair])[0], (pairPoints[pair])[1]
    y1, x1 = (pairPoints[pair])[2], (pairPoints[pair])[3]
    # Compute beta    
    phi1, phi2 = tan(-angleTemplate[y1,x1]), tan(-angleTemplate[y2,x2])
    if 1.0+phi1*phi2 != 0:
        beta = atan((phi1-phi2)/(1.0+phi1*phi2))
    else:
        beta=1.57
    # Compute k
    if x1- refPoint[1] !=0:
        m = atan(-float(y1-refPoint[0])/(x1-refPoint[1]))
    else:
        m =1.57
    k = angleTemplate[y1,x1] - m
    # Scale
    distCentre = sqrt((y1-refPoint[0])*(y1-refPoint[0]) + (x1-refPoint[1])*(x1-refPoint[1]))
    distPoints = sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1))
    c = distCentre/distPoints
    # Insert in the table. The angle is in the interval -pi/2 to pi/2
    entryIndex = int((beta+(pi/2.0))/deltaAngle)
    entry = rTable[entryIndex]
    entry.append((k,c))
            
# Gather evidence for the location
accumulator = createImageF(width, height)
for x1,y1 in itertools.product(range(0, width), range(0, height)):
    if magnitude[y1,x1] != 0:
        # Looking for potential the second points along a line
        secondPoints = [ ]
        m = tan(angle[y1,x1] - alpha)      
        if m>-1 and m<1:
            for delta in range(minimaDistPoints, maxDistPoints):
                x2 = min(x1 + delta, width-1)
                y2 = int(y1 - m * (x2 - x1))
                if y2 > 0 and y2 < height and magnitude[y2,x2] != 0:
                    secondPoints.append((y2,x2))
                    break  
            for delta in range(minimaDistPoints, maxDistPoints):    
                x2 = max(0, x1 - delta)  
                y2 = int(y1 - m * (x2 - x1))
                if y2 > 0 and y2 < height and magnitude[y2,x2] != 0:
                    secondPoints.append((y2,x2))
                    break    
        else:
            m = 1.0/m
            for delta in range(minimaDistPoints, maxDistPoints):
                y2 = min(y1 + delta, height-1)
                x2 = int(x1 - m * (y2 - y1))
                if x2 > 0 and x2 < width and magnitude[y2,x2] != 0:
                    secondPoints.append((y2,x2))
                    break  
            for delta in range(minimaDistPoints, maxDistPoints):
                y2 = max(0, y1 - delta) 
                x2 = int(x1 - m * (y2 - y1))
                if x2 > 0 and x2 < width and magnitude[y2,x2] != 0:
                    secondPoints.append((y2,x2))
                    break  
            
        # Gather evidence  
        numPts = len(secondPoints) > 0
        for ptIndex in range(0, numPts): 
            secondPoint = secondPoints[ptIndex]   
            y2, x2 = secondPoint[0], secondPoint[1]
            distPoints = sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1))

            # Compute beta    
            phi1, phi2 = tan(-angle[y1,x1]), tan(-angle[y2,x2])
            if 1.0+phi1*phi2 != 0:
                beta = atan((phi1-phi2)/(1.0+phi1*phi2))
            else:
                beta=1.57
            # Find entry in table
            entryIndex = int((beta+(pi/2.0))/deltaAngle)
            
            row = rTable[entryIndex]
            numEntriesinRow = len(row)
        
            for kIndex in range(0, numEntriesinRow):
                k, c = (row[kIndex])[0], (row[kIndex])[1]
                distCentre = c * distPoints
                m = tan(angle[y1,x1] - k)                   
                if m>-1 and m<1:
                    for x in range(0, width):
                        y = int(y1 - m * (x - x1))
                        d = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))
                        if y > 0 and y < height and abs(d-distCentre) < 3:
                            accumulator[y,x] += 3.0 - abs(d-distCentre)
                else:
                    m = 1.0/m
                    for y in range(0, height):
                        x = int(x1 - m * (y - y1))
                        d = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))
                        if x > 0 and x < width and abs(d-distCentre) < 3:
                            accumulator[y,x] += 3.0 - abs(d-distCentre)
                    
# Plot accumulator 
plot3DHistogram(accumulator) 

