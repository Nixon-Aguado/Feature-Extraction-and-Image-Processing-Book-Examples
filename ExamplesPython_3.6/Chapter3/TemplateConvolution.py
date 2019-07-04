'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
TemplateConvolution: Filter an image by convolution of a template 
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
imageName = "Giraffe.png"
kernelSize = 5

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create Kernel 
kernelImage = createImageL(kernelSize, kernelSize)

# Set the pixels of a flat kernel
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    kernelImage[y,x] = 1.0  

# Create images to store the result
outputImage = createImageL(width, height)

# Apply kernel
kernelCentre = int((kernelSize - 1) / 2)
for x,y in itertools.product(range(0, width), range(0, height)):
    sumKernel = 0
    sumKernelWeights = 0
    for wx,wy in itertools.product(range(0, kernelSize), range(0, kernelSize)):
        posY = y + wy - kernelCentre
        posX = x + wx - kernelCentre 
        
        if posY > -1 and posY <  height and  posX > -1 and posX <  width:
            sumKernel += inputImage[posY,posX] * kernelImage[wy, wx]
            sumKernelWeights += kernelImage[wy, wx]
    
    if sumKernelWeights > 0:
        outputImage[y,x] = sumKernel / sumKernelWeights
            
# Show output image
showImageL(outputImage)
