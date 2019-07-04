'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
FirstOrderEdgeDetector: Compute gradient by first order derivative  
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
horizEdges = createImageF(width, height)
vertEdges = createImageF(width, height)
outputEdges = createImageF(width, height)

for x,y in itertools.product(range(0, width-1), range(0, height-1)):

    horizEdges[y,x] = abs(float(inputImage[y,x]) - float(inputImage[y+1,x]))
    vertEdges[y,x] = abs(float(inputImage[y,x]) - float(inputImage[y,x+1]))
    
    outputEdges[y,x] = abs(2.0* float(inputImage[y,x]) -                        \
                       float(inputImage[y+1,x]) - float(inputImage[y,x+1]))
                           
# Show horizontal, vertical and all edges
showImageF(horizEdges)
showImageF(vertEdges)
showImageF(outputEdges)
