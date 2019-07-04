'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
HistogramEqualization: Adjust image intensity according to the histogram
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL, createVectorI
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

# Compute histogram of the input image
inputHistogram = computeHistogram(inputImage)

# Vector of integers values to store the number of times a pixel value is repeated
accumulateHistogram = createVectorI(256)

# Create images to store the result
outputImage = createImageL(width, height)

# Distribute the values of the input histogram into the output histogram
sumLevels = 0.0
normalization = float(width * height) / 256
for level in range(0, 256):
    sumLevels += inputHistogram[level]
    accumulateHistogram[level] = sumLevels / normalization
    
# Set the pixels in the output image according to the accumulate histogram
for x,y in itertools.product(range(0, width), range(0, height)):
    outputImage[y,x] = accumulateHistogram[inputImage[y,x]]
        
# Compute histogram of the output image
outputHistogram = computeHistogram(outputImage)
  
# Show output image
showImageL(outputImage)  
    
# Plot histograms
plotHistogram(inputHistogram)
plotHistogram(outputHistogram)
