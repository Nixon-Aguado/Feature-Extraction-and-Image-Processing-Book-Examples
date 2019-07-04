'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 3
FourierConvolution: Filter an image by using the Fourier transform 
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageL, showImageF, createImageF
from FourierUtilities import computeCoefficients, reconstruction, computePowerfromCoefficients
from ImageOperatorsUtilities import imageLogF

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    kernelSize = Size of the kernel
'''
pathToDir = "../../Images/Chapter3/Input/"
imageName = "Eye.png"
kernelSize = 9

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)

# Show input image
showImageL(inputImage)

# Create Kernel 
kernelImage = createImageF(width, height)

# Set the pixels of a flat kernel
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    kernelImage[y, x] = 255.0  
    
# Padding size
widthPad, heightPad = width+kernelSize-1, height+kernelSize-1   

# Padding input
inputPad = createImageF(widthPad, heightPad)
for x,y in itertools.product(range(0, width), range(0, height)):
    inputPad[y,x] = inputImage[y,x] 
 
# Padding and flip template
templatePadFlip = createImageF(widthPad, heightPad)
for x,y in itertools.product(range(0, kernelSize), range(0, kernelSize)):
    templatePadFlip[y, x] = kernelImage[kernelSize-y-1, kernelSize-x-1]
showImageF(templatePadFlip)

# Compute coefficients 
imageCoeff, maxFrequencyW, maxFrequencyH = computeCoefficients(inputPad)
templateCoeff, _, _ = computeCoefficients(templatePadFlip)

# Show the log of the power of the input image and template
powerImage = computePowerfromCoefficients(imageCoeff)
powerImageLog = imageLogF(powerImage)
showImageF(powerImageLog)

powerTemplate = computePowerfromCoefficients(templateCoeff)
powerTemplateLog = imageLogF(powerTemplate)
showImageF(powerTemplateLog)

# Frequency domain multiplication
resultCoeff = createImageF(1 + 2 * maxFrequencyW, 1 + 2 * maxFrequencyH , 2)
for kw,kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),         \
                               range(-maxFrequencyH, maxFrequencyH + 1)):
    w = kw + maxFrequencyW 
    h = kh + maxFrequencyH
    
    resultCoeff[h,w][0] = (imageCoeff[h,w][0] * templateCoeff[h,w][0] -           \
                           imageCoeff[h,w][1] * templateCoeff[h,w][1])
    resultCoeff[h,w][1] = (imageCoeff[h,w][1] * templateCoeff[h,w][0] +           \
                           imageCoeff[h,w][0] * templateCoeff[h,w][1])

# Power result
powerResult = computePowerfromCoefficients(resultCoeff)
powerResultLog = imageLogF(powerResult)
showImageF(powerResultLog)           
                    
# Reconstruction
outputImage = reconstruction(resultCoeff)

outPad = createImageF(width, height)
halfKernel = int(kernelSize/2)
for x,y in itertools.product(range(0, width), range(0, height)):
    outPad[y,x] = outputImage[y + halfKernel, x + halfKernel] 

# Show filter image
showImageF(outPad)
