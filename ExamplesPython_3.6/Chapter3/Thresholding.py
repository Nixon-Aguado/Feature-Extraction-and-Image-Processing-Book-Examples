'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
Thresholding: Create binary image by thresholding 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    threshold = Threshold value
'''
# Data directory
pathToDir = "../../Images/Chapter3/Input/"
imageName = "horse.png"
threshold = 130

# Read image into array
inputImage, width, height  = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create images to store the result
outputImage = createImageL(width, height)

# Set the pixels in the output image
for x,y in itertools.product(range(0, width), range(0, height)):    
    if inputImage[y,x] > threshold:
        outputImage[y,x] = 255
    else:
        outputImage[y,x] = 0
            
# Show output image
showImageL(outputImage)
