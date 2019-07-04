'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
OptimalThresholding: Create binary image by finding an optimal threshold
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createVectorF
from PlotUtilities import plotHistogram
from ImageOperatorsUtilities import computeHistogram, thresholdImage

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

# Create histograms to store cumulative moments
w = createVectorF(256)
m = createVectorF(256)

# Create histograms to store separation 
separability = createVectorF(256)

# Obtain histograms
normalization = 1.0 / float(width * height)
w[0] = normalization * inputHistogram[0]
for level in range(1, 256):
    w[level] = w[level-1] + normalization * inputHistogram[level]
    m[level] = m[level-1] + level * normalization * inputHistogram[level]
    
# Look for the maximum
maximumLevel = 0 
for level in range(0, 256):
    if w[level] * (float(level) - w[level]) != 0:
        separability[level] = float(pow( ( m[255] * w[level] - m[level]), 2)      \
                                   / (w[level] * (float(level) - w[level])))
            
        if separability[level] > separability[maximumLevel]:
            maximumLevel = level

outputImage = thresholdImage(inputImage, maximumLevel)
    
# Show output image
showImageL(outputImage)  
plotHistogram(separability) 
