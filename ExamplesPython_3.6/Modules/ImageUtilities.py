'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
ImageUtilities: Helper module to support image read, create, show and write
                Uses PIL and numpy to work with arrays containing image data  
'''

# Image load by using PIL
from PIL import Image

# Array to store image data
from numpy import array, zeros, amax, amin, clip

# Read an image returning an array with 3 color components
def imageReadRGB(fileName):   
    inputImage = Image.open(fileName)
    width, height = inputImage.size

    inputArray = array(inputImage)
    inputArray = inputArray[:, :, 0:3]

    return inputArray, width, height

# Read an image returning an array with 1 component
def imageReadL(fileName):   
    inputImage = Image.open(fileName)
    width, height = inputImage.size

    inputArray = array(inputImage)
    outputArray = zeros((height, width), dtype='uint8')
    
    # Combine pixel values, for each row and column
    for y in range(0, height):
        for x in range(0, width):  
            rgb = inputArray[y,x]     
            outputArray[y,x] = (int(rgb[0]) + int(rgb[1]) + int(rgb[2])) / 3

    return outputArray, width, height

# Save a image containing float as a gray level image
def imageSaveF(image, fileName, maxScale = 1):   
    height = len(image)
    width = len(image[0])
    
    maximum = amax(image) 
    minimum = amin(image)
          
    # Create a gray level image for display
    outputArray = createImageL(width, height)
    
    # Scale the float values into gray scale values
    for y in range(0, height):
        for x in range(0, width):
            if maximum != minimum:     
                outputArray[y,x] = clip( 255.0 * (image[y,x]*maxScale - minimum) / (maximum - minimum), 0, 255.0)   
    
    # Create output and save
    outputImage = Image.fromarray(outputArray, 'L')
    outputImage.save(fileName)

# Save a image gray level image
def imageSaveL(image, fileName):   
    # Create output and save
    outputImage = Image.fromarray(image, 'L')
    outputImage.save(fileName)
    
# Save a image gray level image
def imageSaveRGB(image, fileName):   
    # Create output and save
    outputImage = Image.fromarray(image, 'RGB')
    outputImage.save(fileName)

# Create a zero array with 3 color components
def createImageRGB(width, height): 
    outputArray = zeros((height, width, 3), dtype='uint8')
    
    return outputArray

# Create a zero array with a color components
def createImageL(width, height): 
    outputArray = zeros((height, width), dtype='uint8')
    
    return outputArray

# Create a zero array two float components
def createImageUV(width, height): 
    outputArray = zeros((height, width, 2), dtype='float')
    
    return outputArray

# Create a zero array two int components
def createImage2I(width, height): 
    outputArray = zeros((height, width, 2), dtype='int')
    
    return outputArray

# Create a zero array with float components
def createImageF(width, height, numComponents = 1): 
    if numComponents == 1:
        outputArray = zeros((height, width), dtype='float')
    else:
        outputArray = zeros((height, width, numComponents), dtype='float')
    
    return outputArray

# Create a multidimensional dimensional array
def createImageNF(*arg): 
    sizes = list(arg)

    # Swap with and height
    if len(arg) > 1:
        
        width = sizes[0]
        sizes[0] = sizes[1]
        sizes[1] = width
    
    outputArray = zeros(sizes, dtype='float')
    
    return outputArray

# Create an array of floats form a list
def createImageFromDataF(dataValues):
    height = len(dataValues)
    width = len(dataValues[0])

    outputArray = zeros((height, width), dtype='float')
    
    # Combine pixel values, for each row and column
    for y in range(0, height):
        for x in range(0, width):     
            outputArray[y,x] = dataValues[y][x]

    return outputArray, width, height

# Copy image scaling gray level
def createScaleImageL(image, grayScale):
    height = len(image)
    width = len(image[0])
    
    outputImage = createImageL(width, height)
    for y in range(0, height):
        for x in range(0, width): 
            outputImage[y,x] = int(grayScale * image[y,x])
    return outputImage

def scaleImageL(image, newWidth, newHeight):
    height = len(image)
    width = len(image[0])

    outputImage= createImageL(newWidth, newHeight)
    
    sW = float(width)/newWidth
    sH = float(height)/newHeight
    
    for y in range(0, newHeight):
        for x in range(0, newWidth):
            xS, yS = int(x * sW), int(y * sH)
            
            outputImage[y,x] = image[yS,xS] 
    
    return outputImage;

# Create a 1D array of floats form a list
def createVectorFromDataF(dataValues):
    width = len(dataValues)
    outputArray = zeros(width, dtype='float')
    
    for x in range(0, width):
        outputArray[x] = dataValues[x]
    
    return outputArray, width

# Create a 1D array of floats 
def createVectorF(width):  
    outputArray = zeros(width, dtype='float')
    
    return outputArray

# Create a 1D array of integers 
def createVectorI(width):
    outputArray = zeros(width, dtype='int')
    
    return outputArray
    
# Show an array containing image data
def showImageRGB(image):
    # Create output and show
    outputImage = Image.fromarray(image, 'RGB')
    outputImage.show()
    
# Show an array containing image data
def showImageL(image):
    # Create output and show
    outputImage = Image.fromarray(image, 'L')
    outputImage.show()

# Show an array containing image data
def showImageF(image, maxScale = 1):
    height = len(image)
    width = len(image[0])
    
    maximum = amax(image) 
    minimum = amin(image)
       
    # Create a gray level image for display
    outputArray = createImageL(width, height)
    
    # Scale the float values into gray scale values
    for y in range(0, height):
        for x in range(0, width):
            if maximum != minimum:  
                outputArray[y,x] = clip( 255.0 * (image[y,x]*maxScale - minimum) / (maximum - minimum), 0, 255.0)  
            else:
                outputArray[y,x] = image[y,x]
    
    # Create output and show
    outputImage = Image.fromarray(outputArray, 'L')
    outputImage.show()