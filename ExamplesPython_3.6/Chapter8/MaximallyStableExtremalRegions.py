'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 8
MaximallyStableRegions: Compute maximally stable regions in an image
'''
        
# Set module functions
from ImageUtilities import imageReadL, showImageL,createImageF, showImageF
from PrintUtilities import printProgress

# Iteration
from timeit import itertools

# Definition to find the most common element in a list
def mostCommon(lst):
    return max(set(lst), key=lst.count)

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    incThreshold = How much the region can increase
    timeThreshold = For how long the region must not grow to be stable
    startL = Start grey level 
    endL = End gray level 
    incL = Gray level increment
    minRegionSize = size of the smallest region
'''
pathToDir = "../../Images/Chapter8/Input/"
imageName = "castle1.png"
incThreshold = 20
timeThreshold = 30
startL = 10
endL = 140
incL = 2
minRegionSize = 50
maxRegionSize = 1000

# Read image into array and show
inputImage, width, height = imageReadL(pathToDir+imageName)
showImageL(inputImage)

# The number of times the region has grown, its size 
timeRegions = { } 
sizeRegions = { }
incSizeRegions = { }

# Regions
regionsImage = createImageF(width, height)

# Stable regions
resultImage = createImageF(width, height)

# Use a threshold to flood regions
nextRegionID = 1
for threshold in range(startL, endL, incL):
    printProgress(threshold - startL,  endL - startL)
    
    # Init the change in size
    for regionID in incSizeRegions:
        incSizeRegions[regionID] = 0
      
    # Repeatedly flood the image to grow regions  
    flooded = True   
    while flooded:
        flooded = False
        growRegion = [ ]
        # For each non-region pixels
        for x,y in itertools.product(range(0, width), range(0, height)):      
            if inputImage[y,x] <= threshold and regionsImage[y,x] == 0:
                
                #  List of neighbours
                n = [ ]
                for wx,wy in itertools.product(range(x-1, x+2), range(y-1, y+2)):
                    if wy>=0 and wy<height and wx>=0 and wx<width:
                        neighbourID = regionsImage[wy, wx]
                        if neighbourID != 0: 
                            n.append(regionsImage[wy, wx])
                
                # Grow the most common
                if(len(n) != 0):
                    mc = mostCommon(n) 
                    growRegion.append((y,x,mc))
                    flooded  = True
        
        for pixel in growRegion:   
            y, x, idRegion = pixel[0] , pixel[1] , pixel[2]               
            regionsImage[y, x] = idRegion
            incSizeRegions[idRegion] += 1
        
    # Repeatedly merge regions
    merged = True
    while merged:
        merged = False
        # For each non-region pixels
        for x,y in itertools.product(range(0, width), range(0, height)):      
            if regionsImage[y,x] != 0:
                
                #  List of neighbours and positions
                n, p = [ ], [ ]
                for wx,wy in itertools.product(range(x-1, x+2), range(y-1, y+2)):
                    if wy>=0 and wy<height and wx>=0 and wx<width:
                        neighbourID = regionsImage[wy, wx]
                        if neighbourID != 0:
                            n.append(regionsImage[wy, wx])
                            p.append((wy, wx))
                
                # Different neighbours, we need to merge
                if len(n) != 0 and len(set(n)) != 1:
                    merged = True
                    
                    unique = set(n)
                    mainRegion = regionsImage[y,x] 
                    if mainRegion in unique:
                        unique.remove(mainRegion)
                    
                    # Merge seeds
                    for otherRegion in p:
                        py, px = otherRegion[0], otherRegion[1]
                        if regionsImage[py, px ] != mainRegion:
                            growRegion.append(( py, px ))
                            regionsImage[py, px ] = mainRegion
                            
                    while len(growRegion) > 0:
                        seed = growRegion.pop()
                        py,px = seed[0], seed[1]
                        
                        regionsImage[py,px] = mainRegion
                        incSizeRegions[mainRegion] += 1    
                        
                        for wx,wy in itertools.product(range(px-1, px+2), range(py-1, py+2)):
                            if wy>=0 and wy<height and wx>=0 and wx<width:
                                if regionsImage[wy,wx] in unique:
                                    regionsImage[wy,wx] = mainRegion
                                    growRegion.append((wy,wx)) 
                    
                    for  regionID in unique:
                        del incSizeRegions[regionID]
                        del timeRegions[regionID]
                        del sizeRegions[regionID]
                        
    # Find new regions
    for x,y in itertools.product(range(0, width), range(0, height)):      
            if inputImage[y,x] <= threshold and regionsImage[y,x] == 0:
                timeRegions[nextRegionID] = 0 
                sizeRegions[nextRegionID] = 0
                incSizeRegions[nextRegionID] = 0

                growRegion = [(y,x)]
                while len(growRegion) > 0:
                    seed = growRegion.pop()
                    py,px = seed[0], seed[1]
                    
                    regionsImage[py,px] = nextRegionID
                    incSizeRegions[nextRegionID] += 1
                    
                    for wx,wy in itertools.product(range(px-1, px+2), range(py-1, py+2)):
                        if wy>=0 and wy<height and wx>=0 and wx<width:
                            if inputImage[wy,wx] <= threshold and regionsImage[wy,wx] == 0:
                                growRegion.append((wy,wx))
                
                nextRegionID += 1
                    
    # Update times for regions
    for idRegion in incSizeRegions:
        # Update the size
        incSize = incSizeRegions[idRegion]
        sizeRegions[idRegion] += incSize
        
        # Update  stable
        if incSize < incThreshold:     
            timeRegions[idRegion] += 1
        else:
            timeRegions[idRegion] = 0
            
        # Stable region condition 
        if timeRegions[idRegion] > timeThreshold and sizeRegions[idRegion] > minRegionSize:
            for x,y in itertools.product(range(0, width), range(0, height)):
                if regionsImage[y,x] == idRegion:
                    resultImage[y,x] = 255 
            timeRegions[idRegion] = 0                
  
showImageF(resultImage)     

                     