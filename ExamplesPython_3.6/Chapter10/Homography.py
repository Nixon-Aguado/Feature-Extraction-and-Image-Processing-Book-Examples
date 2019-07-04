'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 10
Homography: Compute an homography from four corresponding image points and perform the 
            transformation on the image  
'''
        
# Set module functions
from ImageUtilities import imageReadRGB, imageReadL, showImageRGB
from GeometricUtilities import solveSystem, imageTransform
                         
# Math and iteration
from math import sin, cos, sqrt 
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    maskName = Mask image name
'''
pathToDir = "../../Images/Chapter10/Input/"
imageName = "cube1.png"
maskName = "mask1.png"

# Read image data
inputImage, width, height = imageReadRGB(pathToDir + imageName)
maskImage, width, height = imageReadL(pathToDir + maskName)
showImageRGB(inputImage)

# Image centre
centreX, centreY = int(width/2), int(height/2)

# Corresponding points
p = [[116-centreX,202-centreY],[352-centreX,234-centreY],[140-centreX,384-centreY],[344-centreX,422-centreY]]
q = [[118-centreX,168-centreY],[312-centreX,238-centreY],[146-centreX,352-centreY],[322-centreX,422-centreY]]

# Find transform
M = [[-p[0][0], -p[0][1], -1, 0, 0, 0,  p[0][0]*q[0][0], p[0][1]*q[0][0], q[0][0]], \
     [ 0, 0, 0, -p[0][0], -p[0][1], -1, p[0][0]*q[0][1], p[0][1]*q[0][1], q[0][1]], \
     [-p[1][0], -p[1][1], -1, 0, 0, 0,  p[1][0]*q[1][0], p[1][1]*q[1][0], q[1][0]], \
     [ 0, 0, 0, -p[1][0], -p[1][1], -1, p[1][0]*q[1][1], p[1][1]*q[1][1], q[1][1]], \
     [-p[2][0], -p[2][1], -1, 0, 0, 0,  p[2][0]*q[2][0], p[2][1]*q[2][0], q[2][0]], \
     [ 0, 0, 0, -p[2][0], -p[2][1], -1, p[2][0]*q[2][1], p[2][1]*q[2][1], q[2][1]], \
     [-p[3][0], -p[3][1], -1, 0, 0, 0,  p[3][0]*q[3][0], p[3][1]*q[3][0], q[3][0]], \
     [ 0, 0, 0, -p[3][0], -p[3][1], -1, p[3][0]*q[3][1], p[3][1]*q[3][1], q[3][1]], \
     [ 1, 1, 1,     1,        1,     1,        1,               1,            1  ]] 

# Solves the equation A*x=b
b = [0,0,0,0,0,0,0,0,1]
h = solveSystem(M, b)

H = [[h[0], h[1], h[2]],    \
     [h[3], h[4], h[5]],    \
     [h[6], h[7], h[8]] ]
#print(H)

# Transform image and show
tImage = imageTransform(inputImage, maskImage, H)
showImageRGB(tImage)


