'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
PrintUtilities: Helper module to print image data and messages on the standard device  
'''

# Print image range
def printImageRangeRGB(image, rangeWidth, rangeHeight):
    print ("\n")
    for y in range(rangeHeight[0], rangeHeight[1] + 1): 
        print ("[", end=' ')
        for x in range(rangeWidth[0], rangeWidth[1] + 1):
            rgb = image[y,x]
            print ("(", format(rgb[0], '3d'), format(rgb[1], '3d'), format(rgb[2], '3d'), ")", end=' ') 
        print ("]")
        
# Print image range
def printImageRangeL(image, rangeWidth, rangeHeight):
    print ("\n")
    for y in range(rangeHeight[0], rangeHeight[1] + 1): 
        print ("[", end=' ')
        for x in range(rangeWidth[0], rangeWidth[1]  + 1):
            print (format(image[y,x], '3d'), end=' ') 
        print ("]")   
        
# Print image range
def printImageRangeF(image, rangeWidth, rangeHeight, formatValue = ' 3.2f'):
    print ("\n")
    for y in range(rangeHeight[0], rangeHeight[1] + 1): 
        print ("[", end=' ')
        for x in range(rangeWidth[0], rangeWidth[1]  + 1):
            print (format(image[y,x], formatValue), end=' ') 
        print ("]")   
        
# Print text
def printText(text):  
    print (text)
    
# Print text
def printTextSameLine(text):  
    print (text, end=' ')
        
# Some operations can be slow. Print something to show progress
def printProgress(step, totalSteps):
    print (step,"/",totalSteps,"...")
    
