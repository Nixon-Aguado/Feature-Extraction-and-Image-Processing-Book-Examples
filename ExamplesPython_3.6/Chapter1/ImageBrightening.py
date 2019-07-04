'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 1
ImageBrightening: Increase the intensity of an image
'''

# Iteration
from timeit import itertools

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL 
from PrintUtilities import printImageRangeL

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    brightDelta = Increase brightness
    printRange = Image range to print
'''
pathToDir = "../../Images/Chapter1/Input/"
imageName = "Zebra.png"
brightDelta = 80;
printRange = [0, 10]

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Output image
outputImage = createImageL(width, height)

# Set the pixels in the output image
for x,y in itertools.product(range(0, width), range(0, height)):  
    outValue = int(inputImage[y,x]) + brightDelta 
    if outValue < 255:
        outputImage[y,x] = outValue 
    else:
        outputImage[y,x] = 255

# Show images
showImageL(inputImage)
showImageL(outputImage)

# Print image range
printImageRangeL(inputImage, printRange, printRange)
printImageRangeL(outputImage, printRange, printRange)
