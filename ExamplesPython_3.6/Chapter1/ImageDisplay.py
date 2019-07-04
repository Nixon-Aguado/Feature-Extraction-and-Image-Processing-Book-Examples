'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 1
ImageDipslay: Loads and displays an image. Shows a surface and prints pixel data
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF
from PrintUtilities import printImageRangeL
from PlotUtilities import plotColorSurface, plot3DColorHistogram

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter1/Input/"
imageName = "SmoothSquare.png"

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Print pixel's values in an image range
printImageRangeL(inputImage, [0, width-1], [0, height-1])

# Create an image to store the z values for surface
outputZ = createImageF(width, height)

# Three float array to store colors of the surface
colorsRGB = createImageF(width, height, 3)

# Set surface and color values
for x in range(0, width):
    for y in range(0, height):
        pixelValue = float(inputImage[y,x])
        outputZ[y,x] = 255 - pixelValue
        pointColour = float(inputImage[y,x])/255.0
        colorsRGB[y,x] = [pointColour, pointColour, pointColour]

# Plot surface
plotColorSurface(outputZ, colorsRGB, [0,400], 1)

# Plot histogram
plot3DColorHistogram(outputZ, colorsRGB, [0,400])