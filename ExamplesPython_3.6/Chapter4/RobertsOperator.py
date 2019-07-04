'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
RobertsOperator: Compute the image gradient according to  detection by the Roberts operator 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF

# Iteration
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

# Create images to store the results
outputMminus = createImageF(width, height)
outputMplus = createImageF(width, height)
outputEdges = createImageF(width, height)

for x,y in itertools.product(range(0, width-1), range(0, height-1)):

    outputMminus[y,x] = abs(float(inputImage[y,x]) - float(inputImage[y+1,x +1]))
    outputMplus[y,x] = abs(float(inputImage[y+1,x]) - float(inputImage[y,x+1]))  
    outputEdges[y,x] = max(outputMminus[y,x], outputMplus[y,x])
                               
# Show horizontal, vertical and all edges
showImageF(outputMminus)
showImageF(outputMplus)
showImageF(outputEdges)
