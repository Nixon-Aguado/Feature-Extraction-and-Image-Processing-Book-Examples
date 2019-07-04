'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
PrewittOperator: Compute gradient by using the Prewitt operator  
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from PrintUtilities import printImageRangeF
from PlotUtilities import plotQuiver

# Math and iteration
from math import sqrt, atan2
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter4/Input/"
imageName = "Squares.png"

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

outputMagnitude = createImageF(width, height)
outputDirection = createImageF(width, height)

for x,y in itertools.product(range(0, width-1), range(0, height-1)):
    mX,mY = 0.0, 0.0
    for c in range(-1, 2):
        mX += float(inputImage[y - 1, x + c]) - float(inputImage[y + 1, x + c])
        mY += float(inputImage[y + c, x - 1]) - float(inputImage[y + c, x + 1])
    outputMagnitude[y,x] = sqrt(mX * mX + mY * mY)   
    outputDirection[y,x] = atan2(mY, mX)     
            
# Show output image
showImageF(outputMagnitude)
showImageF(outputDirection)

# Print pixel's values in an image range
printImageRangeF(outputDirection, [0, width-1], [0, height-1])

# Plot vectors
plotQuiver(outputMagnitude, outputDirection, 1300)
