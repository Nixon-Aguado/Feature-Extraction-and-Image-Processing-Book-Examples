'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
ImagePropertiesUtilities: Helper module to obtain information from an image
'''

# Array to store image data
from numpy import amax, amin, unravel_index

# Mean form a float image
def meanStddev(image):
    heightImage = len(image)
    widthImage = len(image[0])
    
    m = 0.0
    for x in range(0, widthImage):
        for y in range(0, heightImage):
            m += float(image[y,x])
    m /= float(widthImage) * heightImage
    
    s = 0.0
    for x in range(0, widthImage):
        for y in range(0, heightImage):
            s += (float(image[y,x]) -m) ** 2.0
    
    return m,s

# Return the maximum and minimum
def imageMaxMin(image):   
    maximum = amax(image)
    minimum = amin(image)
    
    return maximum, minimum

# Return the maximum position
def imageArgMax(image):  
    index = unravel_index(image.argmax(), image.shape) 
          
    return index

# Detect a peaks in a 2D image
def peakDetectorImageL(image, peakThreshold, suppWindow = 3):
    peaks = []
    
    height = len(image)
    width = len(image[0])

    for y in range(0, height):
        for x in range(0, width):  
        
            if image[y, x] > peakThreshold:
                peak = True
                for wy in range(y-suppWindow, y+suppWindow+1):
                    for wx in range(x-suppWindow, x+suppWindow+1):
                        if wy>=0 and wy<height and wx>=0 and wx<width  and     \
                            image[y, x] < image[wy, wx]:
                                peak = False
                if peak:
                    peaks.append((y,x)) 
    return peaks

# Detect a peaks in a vector
def peakDetectorVector(image, peakThreshold, suppWindow = 3):
    peaks = []
    
    size = len(image)
    for x in range(0, size):  
        if image[x] > peakThreshold:
            peak = True
            for wx in range(x-suppWindow, x+suppWindow+1):
                if wx>=0 and wx<size and image[x] < image[wx]:
                    peak = False
            if peak:
                peaks.append(x) 
    return peaks



