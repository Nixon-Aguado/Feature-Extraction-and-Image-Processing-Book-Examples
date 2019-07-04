'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 10
Reprojection: Compute a projection from seven corresponding image and 3D points and re-project 
              the image to create a new view of the scene
'''
        
# Set module functions
from ImageUtilities import imageReadRGB, imageReadL, showImageRGB, createImageRGB
from GeometricUtilities import projectionCubePoints, computeProjection, getPointColours, fillImageColours

# Math
from math import sin, cos

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

centreX, centreY = width/2, height/2

# Corresponding points in the cube image and 3D world
pts = [[131-centreX,378-centreY],[110-centreX,188-centreY],
       [200-centreX,70-centreY],[412-centreX,100-centreY],
       [410-centreX,285-centreY],[349-centreX,418-centreY],
       [345-centreX,220-centreY]]

q = [[0,0,1],[0,1,1],
     [0,1,0],[1,1,0],
     [1,0,0],[1,0,1],
     [1,1,1]]

# Obtain the projection
p = computeProjection(pts,q)

# Get the image position of the 3D cube points
npts = 100
xy = projectionCubePoints(npts, p, centreX, centreY)

# Get the colour of the points
colours = getPointColours(xy, maskImage, inputImage)

# Transform the q points and store in qt
qT = [ ]
angY = .3
angX = -.2
for pointNum in range(0,len(q)):
    s = [q[pointNum][0]-.5, q[pointNum][1]-.5, q[pointNum][2]-.5]
    rx = .5 + cos(angY)*s[0] + sin(angY)*s[2]
    ry = .5 + sin(angX)*sin(angY)*s[0] + cos(angX)*s[1] - sin(angX)*cos(angY)*s[2] 
    rz = .5 - cos(angX)*sin(angY)*s[0] + sin(angX)*s[1] + cos(angX)*cos(angY)*s[2] 
    
    qT.append([rx,ry,rz])

# Get the projection of the transformed points
p = computeProjection(pts,qT)

# The position of the cube points according to the projection of the transformed data
xy = projectionCubePoints(npts, p, centreX, centreY)

# Use the colours of the original image and the points of the transformed projection to generate an image
tImage = createImageRGB(width, height)
fillImageColours(colours, xy, tImage)
showImageRGB(tImage)  
