'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
MedianFilter: Noise reduction by median filter 
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
imageName = "Fence.png"
kernelSize = 5

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create images to store the result
outputImage = createImageL(width, height)

# Apply filter
kernelCentre = int((kernelSize - 1) / 2)
for x,y in itertools.product(range(0, width), range(0, height)):
    region = [ ]
    for wx,wy in itertools.product(range(0, kernelSize), range(0, kernelSize)):

        posY = y + wy - kernelCentre
        posX = x + wx - kernelCentre 
        
        if posY > -1 and posY <  height and  posX > -1 and posX <  width:
            region.append(inputImage[posY,posX]) 
    
    numPixels = len(region) 
    if  numPixels > 0:
        region.sort()
        outputImage[y,x] = region[int(numPixels/2)]
        
# Show output image
showImageL(outputImage)
