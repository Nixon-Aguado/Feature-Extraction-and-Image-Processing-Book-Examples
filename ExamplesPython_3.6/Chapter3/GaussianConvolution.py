'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
GaussianConvolution: Filter an image by the convolution of a Gaussian function 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF
from ConvolutionUtilities import applyKernel
from PlotUtilities import plotColorSurface, plot3DColorHistogram
from PrintUtilities import printImageRangeF

# Math and iteration
from math import pow, exp
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    kernelSize = Kernel size
    sigma = Standard deviation
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "Giraffe.png"
kernelSize = 9
sigma = 3.0

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create image to store kernel 
kernelImage = createImageF(kernelSize, kernelSize)

# Three float array to store colors to be used in the surface plot
colorsRGB = createImageF(kernelSize, kernelSize, 3)

# Set the pixels of Gaussian kernel
centre = (kernelSize - 1) / 2
sumValues = 0
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    kernelImage[y,x] = exp( -0.5 * (pow((x - centre)/sigma, 2.0) +    \
                                    pow((y - centre)/sigma, 2.0)) ) 
    sumValues += kernelImage[y,x] 
     
    # This is only set for the plot
    colorsRGB[y,x] = [kernelImage[y,x], kernelImage[y,x], kernelImage[y,x]]

# Normalisation
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    kernelImage[y,x] /= sumValues

# Apply kernel
outputImage =  applyKernel(inputImage, kernelImage)

# Show function
plotColorSurface(kernelImage, colorsRGB)

# Show kernel
plot3DColorHistogram(kernelImage, colorsRGB)

# Print kernel
printImageRangeF(kernelImage, [0, kernelSize-1], [0, kernelSize-1], '7.3')

# Show output image
showImageL(outputImage)
