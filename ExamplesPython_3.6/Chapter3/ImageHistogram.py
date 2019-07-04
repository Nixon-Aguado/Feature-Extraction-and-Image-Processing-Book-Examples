'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
ImageHistogram: Compute the histogram of an image whose gray level values
                are transformed by scale * GrayLevel + Translation 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createVectorI, createImageL
from PlotUtilities import plotHistogram

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    scale = Scale the gray levels
    translation = Value added to the gray levels
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "Horse.png"
scale = 1.0
translation = 0.0

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Transform the image
for x,y in itertools.product(range(0, width), range(0, height)):
    b = int(scale*float(inputImage[y,x]) + translation)
    inputImage[y,x] = max(0, min(b, 255))

# Show transformed image
showImageL(inputImage)

# Vector of integers values to store the number of times a pixel value is repeated
outputHistogram = createVectorI(256)

# Get the number of times a pixel value is found in the image
for x,y in itertools.product(range(0, width), range(0, height)):
    pixelValue = inputImage[y,x]
    outputHistogram[pixelValue] += 1
               
# Plot histogram
plotHistogram(outputHistogram)
