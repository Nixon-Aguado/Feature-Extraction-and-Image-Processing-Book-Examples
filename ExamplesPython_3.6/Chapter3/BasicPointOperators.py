
'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
BasicPointOperators: Applies point operations to an image (sawtooth,logarithmic,exponential) 
                      and show the histogram of the resulting image
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL
from PlotUtilities import plotHistogram
from ImageOperatorsUtilities import computeHistogram

# Math functions and iteration
from math import log, exp
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    intevalSize = Define the sawtooth fixed interval size
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "Horse.png"
intevalSize = 64

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create 3 images to store the result of 3 operators
outputSawtoothImage = createImageL(width, height)
outputLogarithmicImage = createImageL(width, height)
outputExponentialImage = createImageL(width, height)

# Set the pixels in the output image
for x,y in itertools.product(range(0, width), range(0, height)):
    inputValue = int(inputImage[y,x])
    
    # Set the pixels in the sawtooth image
    pixelInInterval = inputValue % intevalSize
    gain = float(pixelInInterval) / float(intevalSize)
    outputSawtoothImage[y,x] = inputValue * gain 
    
    # Set the pixels in the Logarithmic
    outputLogarithmicImage[y,x] = 20 * log(inputValue * 100.0) 
    
    # Set the pixels in the Exponential image
    outputExponentialImage[y,x] = 20 * exp(inputValue / 100.0)  
      
# Compute histograms
histogramSawtoothImage = computeHistogram(outputSawtoothImage)
histogramLogarithmicImage = computeHistogram(outputLogarithmicImage)
histogramExponentialImage = computeHistogram(outputExponentialImage)       
       
# Show output images
showImageL(outputSawtoothImage)
showImageL(outputLogarithmicImage)
showImageL(outputExponentialImage)

# Plot histograms
plotHistogram(histogramSawtoothImage)
plotHistogram(histogramLogarithmicImage)
plotHistogram(histogramExponentialImage)
