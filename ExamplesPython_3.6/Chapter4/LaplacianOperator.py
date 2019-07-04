'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 4
LaplacianOperator: Detect edges by the Laplacian operator
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL, showImageF
from ConvolutionUtilities import createLaplacianKernel, applyKernelF

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    kernelSize = Size of the kernel
    sigma = Standard deviation of the kernel
'''
pathToDir = "../../Images/Chapter4/Input/"
imageName = "Lizard.png"
kernelSize = 12
sigma = 2

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create Kernel 
kernelLaplacian = createLaplacianKernel(kernelSize, sigma)

# Apply kernel
gaussianImage = applyKernelF(inputImage, kernelLaplacian)
    
# Zero-crossing detector  
edges = createImageL(width, height)
kernelCentre = int((kernelSize - 1) / 2)
for x,y in itertools.product(range(1, width-1), range(1, height-1)):
    quadrantValue = [0.0, 0.0, 0.0, 0.0]
    for wx,wy in itertools.product(range(-1, 1), range(-1, 1)):
            quadrantValue[0] += gaussianImage[y+wy, x+wx]
            
    for wx,wy in itertools.product(range(-1, 1), range(0, 2)):
            quadrantValue[1] += gaussianImage[y+wy, x+wx]
    
    for wx,wy in itertools.product(range(0, 2), range(-1, 1)):
            quadrantValue[2] += gaussianImage[y+wy, x+wx]
    
    for wx,wy in itertools.product(range(0, 2), range(0, 2)):
            quadrantValue[3] += gaussianImage[y+wy, x+wx]
                
    maxVal,minVal = max(quadrantValue), min(quadrantValue)      
    
    if maxVal > 0.0 and minVal < 0:
        edges[y,x] = 255  

showImageF(gaussianImage)
showImageL(edges)
