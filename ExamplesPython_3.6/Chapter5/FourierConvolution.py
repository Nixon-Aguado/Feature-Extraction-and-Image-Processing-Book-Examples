'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
Chapter 5
FourierConvolution: Compute the matching of a template in an image by using Fourier convolutions
'''

# Set module functions
from ImageUtilities import imageReadL, showImageL, createImageF, showImageF
from FourierUtilities import  computeCoefficients, reconstruction
from ImagePropertiesUtilities import imageMaxMin
from PrintUtilities import printText
from PlotUtilities import plot3DHistogram 

# Iteration
from timeit import itertools

'''
Parameters:
    pathToDir = Input image directory
    imageName = Input image name
    templateName = Input template image name
    addQuadraticTerm = Set to true to add the square term so that is equivalent to SSD
'''
pathToDir = "../../Images/Chapter5/Input/"
imageName = "Eye.png"
templateName = "EyeTemplate.png"
addQuadraticTerm = True

# Read image into array
inputImage, width, height = imageReadL(pathToDir + imageName)
templateImage, widthTemplate, heightTemplate = imageReadL(pathToDir + templateName)

# We pad the input and template to this size
widthPad = width + widthTemplate - 1
heightPad = height + heightTemplate - 1    

# Pad input 
inputPad = createImageF(widthPad, heightPad)
for x,y in itertools.product(range(0, width), range(0, height)):
    inputPad[y,x] = inputImage[y,x] 

# Pad and invert template
templatePad = createImageF(widthPad, heightPad)
templatePadFlip = createImageF(widthPad, heightPad)
for x,y in itertools.product(range(0, widthTemplate), range(0, heightTemplate)):
    templatePad[y,x] = templateImage[y, x]
    templatePadFlip[y,x] = templateImage[heightTemplate-y-1, widthTemplate-x-1]

# Show input image and template
showImageF(inputPad) 
showImageF(templatePad)  

# Compute correlation in image domain sum of square differences 
squaredTerm = createImageF(widthPad, heightPad)
corrImage = createImageF(widthPad, heightPad)
for x,y in itertools.product(range(0, widthPad), range(0, heightPad)):
    for w,h in itertools.product(range(-widthTemplate+1,1),                         \
                                 range(-heightTemplate+1,1)):
        p, q = x+w, y+h
        if p >=0 and q>=  0 and p < width and q < height:  
            squaredTerm[y,x] += inputPad[q,p] * inputPad[q,p]
            corrImage[y,x] += 2.0 * templatePad[h+heightTemplate-1,w+widthTemplate-1] * inputPad[q,p]
                    
if addQuadraticTerm:   
    for x,y in itertools.product(range(0, widthPad), range(0, heightPad)): 
        corrImage[y,x] += -squaredTerm[y,x]    
    
showImageF(corrImage) 
maxima, minima = imageMaxMin(corrImage) 
plot3DHistogram(corrImage, [2*(minima+maxima)/3, maxima], [15, -47], False)

# Compute Fourier coefficients 
imageCoeff, maxFrequencyW, maxFrequencyH = computeCoefficients(inputPad) 
templateCoeff, _, _ = computeCoefficients(templatePadFlip)

# Frequency domain multiplication defines convolution is space domain
resultCoeff = createImageF(1 + 2 * maxFrequencyW ,1 + 2 * maxFrequencyH , 2)
for kw, kh in itertools.product(range(-maxFrequencyW, maxFrequencyW + 1),         \
                                range(-maxFrequencyH, maxFrequencyH + 1)):
    w = kw + maxFrequencyW 
    h = kh + maxFrequencyH
    resultCoeff[h,w][0] = (imageCoeff[h,w][0] * templateCoeff[h,w][0] -           \
                           imageCoeff[h,w][1] * templateCoeff[h,w][1])
    resultCoeff[h,w][1] = (imageCoeff[h,w][1] * templateCoeff[h,w][0] +           \
                           imageCoeff[h,w][0] * templateCoeff[h,w][1])
                                  
# Inverse Fourier transform  
reconstructedResult = reconstruction(resultCoeff)

# Show convolution 
showImageF(reconstructedResult)
maxima, minima = imageMaxMin(reconstructedResult)
plot3DHistogram(reconstructedResult, [minima, maxima])

# Add square term to define an operator equivalent to SSD
if addQuadraticTerm:
    for x,y in itertools.product(range(0, widthPad), range(0, heightPad)):  
        reconstructedResult[y,x] = -squaredTerm[y,x] + 2.0 * reconstructedResult[y,x]
else:
    for x,y in itertools.product(range(0, widthPad), range(0, heightPad)):  
        reconstructedResult[y,x] = 2.0 * reconstructedResult[y,x]

# Show convolution added the quadratic image term
showImageF(reconstructedResult)  
maxima, minima = imageMaxMin(reconstructedResult) 
plot3DHistogram(reconstructedResult, [2*(minima+maxima)/3, maxima], [15, -47], False)

