'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
Erosion: Erosion morphological filter 
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
imageName = "Logs.png"
kernelSize = 5

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create Kernel 
kernelCentre = (kernelSize - 1) / 2
kernelImage = createImageL(kernelSize, kernelSize)

# Set the pixels of a flat kernel
for x in range(0, kernelSize):
    for y in range(0, kernelSize):
        kernelImage[y,x] =  1

# Create images to store the result
outputImage = createImageL(width, height)

# Apply kernel
kernelCentre = int((kernelSize - 1) / 2)
for x,y in itertools.product(range(0, width), range(0, height)):

    minValue = 256;
    for wx,wy in itertools.product(range(0, kernelSize), range(0, kernelSize)):
        posY = y + wy - kernelCentre
        posX = x + wx - kernelCentre 
        
        if posY > -1 and posY < height and  posX > -1 and posX <  width:
            sub = float(inputImage[posY,posX]) - kernelImage[wy, wx]
           
            if sub > 0 and sub < minValue: minValue = sub 
    
    outputImage[y,x] = minValue
            
# Show output image
showImageL(outputImage)
