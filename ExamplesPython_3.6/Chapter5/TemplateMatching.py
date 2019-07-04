'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 5
TemplateMatching: Compute the matching of a template in an image
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from PrintUtilities import printProgress
from ImagePropertiesUtilities import imageMaxMin
from PlotUtilities import plot3DHistogram 

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    templateName = Input template image name
    thresholdVal = Only pixels in the template with value greater that this are used
                   -1 to use all pixels or 0  to use edges with value >0
'''
pathToDir = "../../Images/Chapter5/Input/"
imageName = "Eye.png"
templateName = "EyeTemplate.png"
thresholdVal = -1 

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)
templateImage, widthTemplate, heightTemplate = imageReadL(pathToDir + templateName)

# Show input image and template
showImageL(inputImage)
showImageL(templateImage)

# Create an accumulator. We look in a reduced size image
accumulator = createImageF(width, height)
  
# Template matching
templateCentreX = int((widthTemplate - 1) / 2)
templateCentreY = int((heightTemplate - 1) / 2)
for x in range(0, width):  
    printProgress(x, width)
    for y in range(0, height):
        for wx,wy in itertools.product(range(0, widthTemplate), range(0, heightTemplate)):
            posY = y + wy - templateCentreY
            posX = x + wx - templateCentreX 
            
            # The threshold is used to accumulate only the edge pixels in an edge template
            # The difference of pixel values is inverted to show the best match as a peak
            if posY > -1 and posY <  height and  posX > -1 and posX <  width and            \
               templateImage[wy,wx] > thresholdVal:
                diff = 1.0 - abs(float(inputImage[posY,posX]) -                             \
                                 float(templateImage[wy, wx])) / 255.0
                
                accumulator[y,x] += diff*diff 
                        
# Show accumulator within a maxima and mininma region
maxima, minima = imageMaxMin(accumulator)
showImageF(accumulator)
plot3DHistogram(accumulator, [minima, maxima])




