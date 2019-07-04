'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
TruncatedMedianFilter: Noise reduction by truncated median filter 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    kernelSize = Size of the kernel
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "artery.png"
kernelSize = 7

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create images to store the result
outputImage = createImageL(width, height)

# Apply filter
kernelCentre = int((kernelSize - 1) / 2)
for x,y in itertools.product(range(0, width), range(0, height)):
        
    # Iterate Window to collect values to compute mean and median 
    region = [ ]
    sumValues = 0
    for wx,wy in itertools.product(range(0, kernelSize), range(0, kernelSize)):
        posY, posX = y + wy - kernelCentre, x + wx - kernelCentre 
        
        if posY > -1 and posY <  height and  posX > -1 and posX <  width:
            sumValues += inputImage[posY,posX]
            region.append(inputImage[posY,posX]) 
    
    # Compute mean and median of the window
    numPixels = len(region)
    if numPixels > 0:
        
        # Mean and median
        mean = sumValues / numPixels 
        region.sort()
        median = region[int(numPixels/2)]
        
        # Upper and low
        upper, lower = 2.0*median-region[0], 2.0*median-region[numPixels-1]
                    
        # Create a list of truncated values 
        truncatedRegion = [ ]            
        for wx,wy in itertools.product(range(0, kernelSize), range(0, kernelSize)):
            posY, posX = y + wy - kernelCentre, x + wx - kernelCentre
            
            if posY > -1 and posY <  height and  posX > -1 and posX <  width:
                if (inputImage[posY,posX] < upper and median < mean) or               \
                   (inputImage[posY,posX] > lower and median > mean):
                    truncatedRegion.append(inputImage[posY,posX])   
                     
        # Compute median of truncated pixels                 
        numTruncatedPixels = len(truncatedRegion) 
        if  numTruncatedPixels > 0:  
            truncatedRegion.sort()
            outputImage[y,x] = truncatedRegion[int(numTruncatedPixels/2)]    
        else:
            outputImage[y,x] = median
    
# Show output image
showImageL(outputImage)