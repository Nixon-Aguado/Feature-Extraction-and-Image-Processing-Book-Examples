'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
ImageRegionUtilities: Helper module to process regions and line segments 
'''

# Images
from ImageUtilities import imageReadL, createImageF, createImageL, showImageL, createVectorF
from ImageOperatorsUtilities import applyCannyEdgeDetector
from ImagePropertiesUtilities import imageMaxMin

# Math and iteration
from math import pi, sqrt, sin, cos, atan, atan2, factorial, exp
from random import shuffle, sample 
from timeit import itertools

# Array to store image data
from numpy import amax, amin
import numpy as np


# Combination function
def nCr(n, r):
    if n < 0 or r < 0 or n-r < 0: return 0
    return factorial(n) / (factorial(r) * factorial(n-r))

# Rising factorial function
def risingFactorial(x, n):
    p = 1
    for k in range(0,n):
        p = p * (x + k) 
    return p
   
# Compute reference point in an edge image
def computeReferencePoint(edgeImage):
    height, width = len(edgeImage), len(edgeImage[0]) 

    refPoint = [0,0]
    edgePoints = []
    for x,y in itertools.product(range(0, width), range(0, height)):
        if edgeImage[y,x] != 0:
            refPoint[0] += y
            refPoint[1] += x
            edgePoints.append((y,x))
    numPts = len(edgePoints)
    refPoint = [int(refPoint[0]/numPts),int(refPoint[1]/numPts)]
    
    return refPoint, edgePoints

# Return the longest segment from edges
def findLongestSegment(edges):
    height, width = len(edges), len(edges[0])
    # Find line segments 
    segmentsList = []
    segmentsImage = createImageF(width, height)
    maxSegmentLenght = 0
    maxSegmentIndex = 0
    for x,y in itertools.product(range(0, width), range(0, height)):
        if edges[y,x] != 0 and segmentsImage[y,x] == 0:
            segment = [ ]
            segmentPoints = [(y,x)]
            segmentsImage[y,x] = 255
            while len(segmentPoints) > 0:
                yc = (segmentPoints[0])[0]
                xc = (segmentPoints[0])[1]
                segment.append((yc,xc))
                segmentPoints = segmentPoints[1:]
                
                for dx,dy in itertools.product(range(-1,2), range(-1,2)):
                    xn, yn = xc+dx, yc+dy
                    if dx!=0 or dy!=0 and xn > 0 and yn > 0 and xn < width and yn < height:
                        if edges[yn,xn] != 0 and segmentsImage[yn,xn] == 0:
                            segmentPoints.append((yn,xn))
                            segmentsImage[yn,xn] = 255
        
            segmentsList.append(segment)
            if len(segment) > maxSegmentLenght:
                maxSegmentLenght = len(segment)   
                maxSegmentIndex = len(segmentsList) - 1 
                  
    mainSegment = [] 
    segment = segmentsList[maxSegmentIndex]
    curentElement = segment.pop(0)
    sy,sx = curentElement[0], curentElement[1]

    mainSegment.append(curentElement)

    numPoints = len(segment)
    while numPoints > 0:
        closestElement = [0, float("inf")]
        cy,cx = curentElement[0], curentElement[1]
        for p in range(0, numPoints):
            y,x = (segment[p])[0], (segment[p])[1]
            d = sqrt((cx-x) * (cx-x)  + (cy-y)  * (cy-y) ) 
            if d < closestElement[1] or (d == closestElement[1] and y > cy):
                closestElement = [p, d]       
        
        # If we are closer to the first point, then end now
        dFirst = sqrt((cx-sx) * (cx-sx)  + (cy-sy)  * (cy-sy) ) 
        if (cx!=sx or cy!=sy) and 2*dFirst < closestElement[1]:
            break
            
        curentElement = segment.pop(closestElement[0])
        numPoints = len(segment)
    
        mainSegment.append(curentElement)
    numPoints = len(mainSegment) 
    
    # Average to get more accurate direction
    averageSize = 1
    totalPixels = float(1 + 2*averageSize)
    mainSegmentAverage = [ ]
    for p in range(0, numPoints): 
        y,x = 0, 0
        for w in range(-averageSize, averageSize+1):
            p1 = p + w
            if p1 < 0: p1 = p1 + numPoints
            if p1 >= numPoints: p1 = p1 - numPoints
            x += (mainSegment[p1])[1]
            y += (mainSegment[p1])[0]
        mainSegmentAverage.append((y/totalPixels, x/totalPixels))
    
    return mainSegmentAverage

# Return the longest segment from an image
def findLongestCentredSegmentinImage(imageName, gaussianKernelSize, sobelKernelSize, upperT, lowerT):
    # Read image into array and show
    inputImage, width, height = imageReadL(imageName)

    # Compute edges and find the segment in the image
    magnitude, _ = applyCannyEdgeDetector(inputImage, gaussianKernelSize, sobelKernelSize, upperT, lowerT)
    mainSegmentAverage = findLongestSegment(magnitude)
    
    # Compute centre
    numPoints = len(mainSegmentAverage)
    centre = [0,0]
    for p in range(0, numPoints):
        centre[0] += (mainSegmentAverage[p])[0]
        centre[1] += (mainSegmentAverage[p])[1]
    centre[0] /= numPoints
    centre[1] /= numPoints    
      
    # Respect to the center and convert to an image array
    shape= createImageF(numPoints, 2)
    for p in range(0, numPoints):
        y,x = (mainSegmentAverage[p])[0], (mainSegmentAverage[p])[1]
        shape[0, p] = y-centre[0]
        shape[1, p] = x-centre[1]
    
    return centre, shape, width, height

def findLongesSegmentinImage(imageName, gaussianKernelSize, sobelKernelSize, upperT, lowerT):
    # Read image into array and show
    inputImage, width, height = imageReadL(imageName)

    # Compute edges and find the segment in the image
    magnitude, _ = applyCannyEdgeDetector(inputImage, gaussianKernelSize, sobelKernelSize, upperT, lowerT)
    mainSegmentAverage = findLongestSegment(magnitude)
      
    # Convert to an image array
    numPoints = len(mainSegmentAverage)
    shape= createImageF(numPoints, 2)
    for p in range(0, numPoints):
        y,x = (mainSegmentAverage[p])[0], (mainSegmentAverage[p])[1]
        shape[0, p] = y
        shape[1, p] = x
    
    return shape, width, height

# Get a list with the pixels outside a backgroundRange
def pixlesList(image, backgroundRange):
    listPixels = [ ]
    height, width = len(image), len(image[0])
    for x,y in itertools.product(range(0, width), range(0, height)):
        if image[y,x] < backgroundRange[0] or image[y,x]  > backgroundRange[1]: 
            listPixels.append((y,x,1))
    return listPixels

def edgesList(image, shapeImage, backgroundRange):
    edgePixels = [ ]
    height, width = len(image), len(image[0])
    numPoints = len(shapeImage)
    for indexPixel in range(0, numPoints):
        y, x = (shapeImage[indexPixel])[0], (shapeImage[indexPixel])[1]
        edge = False
        for wx,wy in itertools.product(range(-1, 2), range(-1, 2)):        
            posX, posY = x + wx, y+ wy
            if posY > -1 and posY < height and  posX > -1 and posX < width:
                if image[posY,posX] >= backgroundRange[0] and image[posY,posX] <= backgroundRange[1] :
                    edge = True
        if edge:
            edgePixels.append((y,x))
    return edgePixels


def computeAngularFunctions(shape):
    # Compute the accumulative arc lengths 
    numPoints = len(shape[0])
    sumArcLenghts = []  
    y0, x0 = shape[0, numPoints-1], shape[1, numPoints-1]
    shapeLenght = 0.0
    for p in range(0, numPoints):
        y,x = shape[0,p], shape[1,p]
        shapeLenght += sqrt((y-y0)*(y-y0) + (x-x0)*(x-x0))
        sumArcLenghts.append(shapeLenght)
        y0,x0 = y,x
  
    # Normalized lengths
    normArcLenghts = []
    for p in range(0, numPoints):
        normArcLenghts.append((2.0*pi*sumArcLenghts[p])/shapeLenght);

    # Compute angular function by an average window
    windowSize = [5,10]
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
        diff = angle-angle0
        if diff < pi:
            diff += 2.0* pi
        if diff > pi:
            diff -= 2.0 * pi
        sumAngle += diff
        cumulativeFunc.append(sumAngle)
        angle0 = angle    
    
    # Compute cumulative angular accumulated
    cumulativeNormFunc = [ ]
    for p in range(0, numPoints):
        cumulativeNormFunc.append(cumulativeFunc[p]+normArcLenghts[p])
    
    return sumArcLenghts, normArcLenghts, angularFunc, cumulativeFunc, cumulativeNormFunc

def weightedKrawtchoukPolynomials(p, width):
    
    # Data containers
    sigma = createVectorF(width)
    ro = createVectorF(width)
    K = createImageF(width,width)
    
    # Coefficient size
    N = width-1

    # Weight
    for x in range(0,width): 
        sigma[x] = nCr(N, x) * pow(p,x) * pow(1-p,N-x)
        
    # Scale factor. Commented direct computation and using for to avoid factorial
    #for n in range(0,width): 
    #    ro[n] = pow(-1,n) * pow((1-p)/p,n) * (float(factorial(n)) / risingFactorial(-N, n))
    ro[0] = 1
    for n in range(1,N): 
        ro[n] = (-1*((1.0-p)/p)*n/(-N+(n-1)))*ro[n-1]
    ro[N]=(((1.0-p)/p)*N)*ro[N-1]
          
    # Krawtchouk matrix that store result of the polynomial
    # Each row is a polynomial each column is the polynomial value for an x value
    # Alternatively, we could have used the polynomial generating function
    q = 1.0/p
    for n,x in itertools.product(range(0, width), range(0, width)):
        for s in range(0,width):
            K[n,x] += pow(-1,s) * nCr(N-x, n-s) * nCr(x, s) * pow(q-1,n-s) 
                
    # Normalize rows for stability 
    for n in range(0,width):     
        scale = K[n,0]
        for x in range(0,width):
            K[n,x] /= scale          
    
    # Obtain the coefficients A of the polynomials from K
    # Solve for the coefficients A in A*C = K
    C = createImageF(width,width)
    for n,x in itertools.product(range(0, width), range(0, width)):
        C[n,x] = pow(x,n)
        
    CT = np.transpose(C)
    KT = np.transpose(K)
    AT = np.linalg.solve(CT, KT) # solves the equation A*x=b   A*C = k, C'*A' = K'
    A = np.transpose(AT)

    # Product defining the weighted
    w = createImageF(width,width)
    for n,x in itertools.product(range(0, width), range(0, width)):
        w[n,x] = sqrt(sigma[x]/ro[n])
        
    return K, A, sigma, ro, w

def geometricMoments(pixelList, numMoments):
    numPoints = len(pixelList)
      
    # Compute moments
    M = createImageF(numMoments,numMoments)
    for m,n in itertools.product(range(0, numMoments), range(0, numMoments)):
        for indexPixel in range(0, numPoints):
            y = (pixelList[indexPixel])[0]
            x = (pixelList[indexPixel])[1]
            val = (pixelList[indexPixel])[2]
            M[n,m] += (x**n) * (y**m) * val
            
    return M

def geometricInvariantMoments(pixelList, numMoments):
    numPoints = len(pixelList)
      
    # Compute moments
    M = createImageF(numMoments,numMoments)
    for m,n in itertools.product(range(0, numMoments), range(0, numMoments)):
        for indexPixel in range(0, numPoints):
            y = (pixelList[indexPixel])[0]
            x = (pixelList[indexPixel])[1]
            val = (pixelList[indexPixel])[2]
            M[n,m] += (x**n) * (y**m) * val
            
    # Geometric central Moments
    xc,yc = M[1,0]/M[0,0], M[0,1]/M[0,0]
    
    m11 = M[1,1]/M[0,0] - xc*yc
    m20 = M[2,0]/M[0,0] - xc**2
    m02 = M[0,2]/M[0,0] - yc**2
    
    if m20 < m02:
        t = -(0.5 * atan(2.0*m11/(m20-m02)) + pi/2.0)
    else:
        t = -(0.5 * atan(2.0*m11/(m20-m02)))
    
    # Geometric invariant moments
    v = createImageF(numMoments,numMoments)
    vn = createImageF(numMoments,numMoments)
    for m,n in itertools.product(range(0, numMoments), range(0, numMoments)):
        for indexPixel in range(0, numPoints):
            y = (pixelList[indexPixel])[0]
            x = (pixelList[indexPixel])[1]
            val = (pixelList[indexPixel])[2]
            v[n,m] += ((x-xc)*cos(t) - (y-yc)*sin(t))**n * ((x-xc)*sin(t) + (y-yc)*cos(t))**m  * val
        l = (1 + ((n + m) / 2.0))
        vn[n,m] = v[n,m] / pow(M[0,0],l)
                
    return vn


# WaterShed transform
def watherShed(distanceImage, shapeImage, suppWindow):
    height, width = len(distanceImage), len(distanceImage[0])
    watershedImage = createImageF(width, height)
    
        # Initial regions by finding the maximum
    regionIndex = 1   # Start id for a region. Any number different from zero
    numPoints = len(shapeImage)
    for indexPixel in range(0, numPoints):        
        y, x = (shapeImage[indexPixel])[0], (shapeImage[indexPixel])[1]
        if watershedImage[y,x] == 0:
            peak = True
            for wx,wy in itertools.product(range(x-suppWindow, x+suppWindow+1), \
                                           range(y-suppWindow, y+suppWindow+1)):
                if wy>=0 and wy<height and wx>=0 and wx<width:               
                    if watershedImage[wy, wx] != 0 or                           \
                       distanceImage[y, x] < distanceImage[wy, wx]:
                        peak = False
            if peak:
                for wx,wy in itertools.product(range(x-suppWindow, x+suppWindow+1), \
                                               range(y-suppWindow, y+suppWindow+1)):
                    if wy>=0 and wy<height and wx>=0 and wx<width:
                        watershedImage[wy, wx] = regionIndex
                regionIndex += 1
                
    floodRegion = [ ] # The region we need to flood
    for indexPixel in range(0, numPoints):        
        y, x = (shapeImage[indexPixel])[0], (shapeImage[indexPixel])[1]
        if watershedImage[y,x] == 0:
            floodRegion.append((y,x))
                
    # This is not required. We do it to get a better display
    # Create random regions ID. We change the ID for a random value so we get a random gray level when showing the regions 
    c = sample(range(regionIndex), regionIndex) 
    for indexPixel in range(0, numPoints):        
        y, x = (shapeImage[indexPixel])[0], (shapeImage[indexPixel])[1]
        if watershedImage[y, x] != 0:
            watershedImage[y, x] = c[int(watershedImage[y, x])] + 1
    
    # Flooding
    maxDistance, _ = imageMaxMin(distanceImage)
    for floodValue in range(int(maxDistance), 0, -1):
        flooded = True
        while flooded:
            flooded = False
            newFloodRegion = [ ]
            growRegion = [ ]
            shuffle(floodRegion)
            for indexPixel in range(0, len(floodRegion)):
                y, x = (floodRegion[indexPixel])[0], (floodRegion[indexPixel])[1]
                
                # Points not flooded will be considered in following iterations
                if distanceImage[y,x] <= floodValue:
                    newFloodRegion.append((y,x))
                else:
                    #  list of neighbours
                    n = [ ]
                    for wx,wy in itertools.product(range(-1, 2), range(-1, 2)):        
                        posX, posY = x + wx, y+ wy
                        if posY > -1 and posY < height and  posX > -1 and posX < width:
                            if watershedImage[posY, posX] != 0:
                                n.append(watershedImage[posY, posX])
                     
                    # No neighbours, so we cannot grow
                    if(len(n) == 0):
                        newFloodRegion.append((y,x))
                    else:
                        # Grow of only one type of region
                        if len(set(n)) == 1:
                            growRegion.append((y,x,n[0]))
                            flooded  = True
                            
            for pixel in growRegion:                   
                watershedImage[pixel[0], pixel[1]] = pixel[2] 
            
            floodRegion = newFloodRegion
      
    # Set the borders
    shedID = regionIndex + 1
    for indexPixel in range(0, numPoints):        
        y, x = (shapeImage[indexPixel])[0], (shapeImage[indexPixel])[1]
        if watershedImage[y,x] == 0 and distanceImage[y, x] > 0.5:
            watershedImage[y, x] = shedID    
               
    return watershedImage


# Return the maximum and minimum points in a shape
def shapeMaxMin(shape):  
    x,y = shape[1,:], shape[0,:]  
      
    maxY, maxX = amax(x), amax(y)
    if maxY > maxX:
        maxX = maxY
    
    minY = amin(x)
    minX = amin(y)
    if minY < minX:
        minX = minY
        
    return maxX, minX

def showShapeinImage(shape, centre, width, height):
    segmentsImage = createImageL(width, height)
    numPoints = len(shape[0])
    for p in range(0, numPoints):
        y,x = int(centre[0]+shape[0,p]), int(centre[1]+shape[1,p])
        if x > 0 and y > 0 and x < width and y < height:
            segmentsImage[y,x] = 255
    showImageL(segmentsImage)
    
 
# Compute density function from a region in an image
def densityHistogram(image, position, regionRadius, sigma, histoSize):
    height = len(image)
    width = len(image[0])  
    
    # Quantization scale
    colourScale = 256.0 / histoSize

    histogram = createImageF(histoSize, histoSize)
    sumValue = 0
    for deltaX, deltaY in itertools.product(range(-regionRadius[0],regionRadius[0]), range(-regionRadius[1], regionRadius[1])):
    
        x, y  = position[0] + deltaX, position[1] + deltaY    
        if x>0 and y>0 and x<width and y<height :
    
            w = exp(-(deltaX*deltaX + deltaY*deltaY)/(2*sigma*sigma));      
            rgb = image[y,x] / 256.0  
  
            Cb = int((128 - 37.79*rgb[0] - 74.203*rgb[1] +    112*rgb[2])/colourScale)
            Cr = int((128 +   112*rgb[0] - 93.786*rgb[1] - 18.214*rgb[2])/colourScale)   
        
            histogram[Cr,Cb] += w
            sumValue += w

    for r,b in itertools.product(range(0, histoSize), range(0, histoSize)):
        histogram[r,b] /= sumValue

    return histogram

# Get a 2D colour description from a RGB value
def colourFeature(rgb, colourScale):
    nRGB = rgb / 256.0  
    cB = int((128 - 37.79*nRGB[0] - 74.203*nRGB[1] +    112*nRGB[2])/colourScale)
    cR = int((128 +   112*nRGB[0] - 93.786*nRGB[1] - 18.214*nRGB[2])/colourScale) 
    
    return cB, cR

# Implementation of meanshift
def meanShift(inputImage, q, sizeReg, sigma, histoSize, newPos):
    # Weights
    weights = createImageF(2*sizeReg[0], 2*sizeReg[1])
    
    currPos = [0, 0]
    
    colourScale = 256.0 / histoSize

    while(currPos != newPos):
        currPos = newPos
        qs = densityHistogram(inputImage, currPos, sizeReg, sigma, histoSize)

        # Weights
        for deltaX, deltaY in itertools.product(range(-sizeReg[0],sizeReg[0]),   \
                                                range(-sizeReg[1], sizeReg[1])):
                    
            # Position of the pixel in the image and in the weight array
            x, y = currPos[0] + deltaX, currPos[1] + deltaY
            px,py = deltaX+sizeReg[0], deltaY+sizeReg[1] 
            
            # Features
            Cb,Cr= colourFeature(inputImage[y,x], colourScale)
              
            # Update
            if qs[Cr, Cb] > 0:
                weights[py, px] = sqrt(q[Cr, Cb] / qs[Cr, Cb])
            else:
                weights[py, px] = sqrt(q[Cr, Cb] / .000000000001)
    
        # Compute mean shift sums
        meanSum = [0, 0]
        kernelSum = 0
        for deltaX, deltaY in itertools.product(range(-sizeReg[0],sizeReg[0]),   \
                                                range(-sizeReg[1], sizeReg[1])):

            # Position of the pixel in the image
            x, y  = currPos[0] + deltaX, currPos[1] + deltaY 

            # Kernel parameter 
            w = exp(-(deltaX*deltaX + deltaY*deltaY)/(2*sigma*sigma));      

            # Weight index
            px, py = deltaX+sizeReg[0], deltaY+sizeReg[1]
   
            # Mean sum
            meanSum[0] +=  w * weights[py, px] * x
            meanSum[1] +=  w * weights[py, px] * y
        
            # Kernel sum
            kernelSum += w * weights[py, px]
        
        # Mean shift 
        newPos = [int(meanSum[0] / kernelSum), int(meanSum[1] / kernelSum)]
  
    return newPos

# Back project a source image into the target
def backProjection(sourceImage, targetImage, qSource, pSource, pTarget, sizeReg, histoSize):
    height, width = len(sourceImage), len(sourceImage[0])  
    
    colourScale = 256.0 / histoSize
    
    # Projection
    projectionSource = createImageF(width, height)
    projectionTarget = createImageF(width, height)
    for x, y in itertools.product(range(0,width), range(0, height)):
        
        Cb,Cr = colourFeature(sourceImage[y,x], colourScale)  
        projectionSource[y,x] = qSource[Cr,Cb]
           
        Cb,Cr = colourFeature(targetImage[y,x], colourScale)    
        projectionTarget[y,x] = qSource[Cr,Cb]
            
    # Compute geometric moments
    momS = createImageF(3, 3)
    momT = createImageF(3, 3)
    sizeSearch = [int(sizeReg[0] *1.5), int(sizeReg[1] *1.5)]
    for deltaX, deltaY in itertools.product(range(-sizeSearch[0], sizeSearch[0]),   \
                                            range(-sizeSearch[1], sizeSearch[1])):
        x, y  = pSource[0] + deltaX, pSource[1] + deltaY 
        for m,n in itertools.product(range(0, 3), range(0, 3)):
            momS[n,m] += (x**n) * (y**m) * projectionSource[y,x]
            
        x, y  = pTarget[0] + deltaX, pTarget[1] + deltaY 
        for m,n in itertools.product(range(0, 3), range(0, 3)):
            momT[n,m] += (x**n) * (y**m) * projectionTarget[y,x]    
            
    xc,yc = momS[1,0]/momS[0,0], momS[0,1]/momS[0,0]
    a = momS[2,0]/momS[0,0] - xc*xc;
    b = 2*(momS[1,1]/momS[0,0] - xc * yc);
    c = momS[0,2]/momS[0,0]- yc*yc;
    sxS = int(sqrt((a+c-sqrt(b*b+(a-c)*(a-c))/2))); 
    syS = int(sqrt((a+c+sqrt(b*b+(a-c)*(a-c))/2)));       
            
    xc,yc = momT[1,0]/momT[0,0], momT[0,1]/momT[0,0]
    a = momT[2,0]/momT[0,0] - xc*xc;
    b = 2*(momT[1,1]/momT[0,0] - xc * yc);
    c = momT[0,2]/momT[0,0]- yc*yc;
    sx = int(sqrt((a+c-sqrt(b*b+(a-c)*(a-c))/2)));
    sy = int(sqrt((a+c+sqrt(b*b+(a-c)*(a-c))/2)));
    
    sy = sy * sizeReg[1] /  syS
    sx = sx * sizeReg[0] /  sxS
    
    return [int(xc),int(yc)], [int(sx),int(sy)]

# Back project an image
def backProjectionImage(image, q, histoSize):    
    height, width = len(image), len(image[0])  
    colourScale = 256.0 / histoSize
    
    # Projection
    projection = createImageF(width, height)
    for x, y in itertools.product(range(0,width), range(0, height)):
        
        Cb,Cr = colourFeature(image[y,x], colourScale)  
        projection[y,x] = q[Cr,Cb]
         
    return projection

# Determine the size of a region
def regionSize(backProjImage, newBackProjImage, pos, newPos, sizeReg):
    # Compute geometric moments
    momS = createImageF(3, 3)
    momT = createImageF(3, 3)
    sizeSearch = [int(sizeReg[0] *1.5), int(sizeReg[1] *1.5)]
    for deltaX, deltaY in itertools.product(range(-sizeSearch[0], sizeSearch[0]),   \
                                            range(-sizeSearch[1], sizeSearch[1])):
        x, y  = pos[0] + deltaX, pos[1] + deltaY 
        for m,n in itertools.product(range(0, 3), range(0, 3)):
            momS[n,m] += (x**n) * (y**m) * backProjImage[y,x]
            
        x, y  = newPos[0] + deltaX, newPos[1] + deltaY 
        for m,n in itertools.product(range(0, 3), range(0, 3)):
            momT[n,m] += (x**n) * (y**m) * newBackProjImage[y,x]    
            
    xc,yc = momS[1,0]/momS[0,0], momS[0,1]/momS[0,0]
    a = momS[2,0]/momS[0,0] - xc*xc;
    b = 2*(momS[1,1]/momS[0,0] - xc * yc);
    c = momS[0,2]/momS[0,0]- yc*yc;
    sxS = int(sqrt((a+c-sqrt(b*b+(a-c)*(a-c))/2))); 
    syS = int(sqrt((a+c+sqrt(b*b+(a-c)*(a-c))/2)));
           
    xc,yc = momT[1,0]/momT[0,0], momT[0,1]/momT[0,0]
    a = momT[2,0]/momT[0,0] - xc*xc;
    b = 2*(momT[1,1]/momT[0,0] - xc * yc);
    c = momT[0,2]/momT[0,0]- yc*yc;
    sx = int(sqrt((a+c-sqrt(b*b+(a-c)*(a-c))/2)));
    sy = int(sqrt((a+c+sqrt(b*b+(a-c)*(a-c))/2)));
    
    sy = sy * sizeReg[1] /  syS
    sx = sx * sizeReg[0] /  sxS
    
    return [int(xc),int(yc)], [int(sx),int(sy)]
