'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
IntensityNormalization: Normalize an image
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL
from ImagePropertiesUtilities import imageMaxMin
from PlotUtilities import plotHistogram
from ImageOperatorsUtilities import computeHistogram

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "Horse.png"

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create image to store the normalization
outputNormalizedImage = createImageL(width, height)

# Maximum and range
maxVal, miniVal = imageMaxMin(inputImage) 
brightRange = float(maxVal - miniVal) 

# Set the pixels in the output image
for x,y in itertools.product(range(0, width), range(0, height)):
    
    # Normalize the pixel value according to the range
    outputNormalizedImage[y,x] = round((inputImage[y,x] - miniVal) * 255.0 / brightRange)
      
# Compute histogram
histogramNormalizedImage = computeHistogram(outputNormalizedImage)
           
# Show output image and plot histogram
showImageL(outputNormalizedImage)
plotHistogram(histogramNormalizedImage)
