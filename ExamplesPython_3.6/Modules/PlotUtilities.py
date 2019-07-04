'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
PlotUtilities: Helper module to show plots of image data   
          Uses numpy  and matplotlib to work with arrays and plots  
'''

# Math functions
import math

# Array to store data
import numpy as np

# Array to store image data
from numpy import zeros

# Plot functions
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Two alternative ways to import Axes3D.  
import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D
#from mpl_toolkits.mplot3d import Axes3D

# Plot a histogram of data
def plotHistogram(data, plotRange = [0, 0], barSepartion = 1):
     
    width = len(data)
    x = np.linspace(0, width-1, width)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca()
    axes.bar(x, data, barSepartion)
    
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    
    if plotRange[0] != 0 or plotRange[1] != 0:
        axes.set_ylim([plotRange[0], plotRange[1]])

    plt.show()
    
# Plot a histogram of data
def plotCurve(data, rangeY = [0, 0], rangeX = [0, 0]):
    
    width = len(data)
    x = np.linspace(0, width-1, width)
    
    if rangeX[0] != 0 or rangeX[1] != 0:
        x = np.linspace(rangeX[0], rangeX[1], width)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca()
    axes.plot(x, data)
    
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    
    if rangeY[0] != 0 or rangeY[1] != 0:
        axes.set_ylim([rangeY[0], rangeY[1]])


    plt.show()
    
def plot2Curves(data1, data2, rangeY = [0, 0]):
    
    width = len(data1)
    x = np.linspace(0, width-1, width)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca()
    axes.plot(x, data1, marker='o' )
    axes.plot(x, data2, marker='o' )
    
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    
    if rangeY[0] != 0 or rangeY[1] != 0:
        axes.set_ylim([rangeY[0], rangeY[1]])
        
    plt.show()
    
def plotCurves(data, rangeY = [0, 0]):
    
    width = data.shape[1]
    height = data.shape[0]
    
    x = np.linspace(0, width-1, width)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca()
    
    for curveNum in range(1, height):
        axes.plot(x, data[curveNum,:])
    
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    
    if rangeY[0] != 0 or rangeY[1] != 0:
        axes.set_ylim([rangeY[0], rangeY[1]])

    plt.show()
  
def plotCurveXY(dataX, dataY, rangeY = [0, 0]):
    if rangeY[0] == 0 and rangeY[1] == 0:
        rangeY = [min(dataY), max(dataY)]
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca()
    axes.plot(dataX, dataY)
    
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    
    deltaRange = (rangeY[1] - rangeY[0])/20.0 
    if rangeY[0] != 0 or rangeY[1] != 0:
        axes.set_ylim([rangeY[0] -deltaRange, rangeY[1]+deltaRange])

    plt.show()

# Plot a histogram of data
def plot3DColorHistogram(dataZ, colorsRGB, zRange = [0, 0]): 
    
    width = dataZ.shape[1]
    height = dataZ.shape[0]
                         
    # Convert data to array                     
    zData = np.array(dataZ)
    
    # Create an X-Y mesh of the same dimension as the 2D data
    xData, yData = np.meshgrid( np.arange(width), np.arange(height))
    
    # Flatten the arrays so that they may be passed to "axes.bar3d".
    xData = xData.flatten()
    yData = yData.flatten()
    zData = zData.flatten()
    
    colorsArray = colorsRGB.reshape((width * height, 3)) 
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca(projection='3d')
    axes.bar3d( xData, yData, np.zeros(len(zData)), .98, .98, zData, color=colorsArray, alpha=1.0, zsort='max') 
    
    axes.set_xlim3d(0, width)
    axes.set_ylim3d(0, height)
    if zRange[0] != 0 or zRange[1] != 0:
        axes.set_zlim3d([zRange[0], zRange[1]])

    plt.show()
    
# Plot a histogram of data
def plot3DHistogram(dataZ, zRange=[0, 0], axesView=[0,0], zTicks = True): 
    
    # Function formatter to simulate zlim_3d
    def major_formatter(x, pos):
        return "{:.1f}".format(x+zRange[0])
    
    width = dataZ.shape[1]
    height = dataZ.shape[0]
    
    # Z clipping fails to clip bar charts
    if zRange[0] != 0:
        for x in range(0, width):
            for y in range(0, height):
                dataZ[y,x] -= zRange[0]
                if dataZ[y,x] < 0:
                    dataZ[y,x] = 0
                         
    # Convert data to array                     
    zData = np.array(dataZ)
    
    # Create an X-Y mesh of the same dimension as the 2D data
    xData, yData = np.meshgrid( np.arange(width), np.arange(height))
    
    # Flatten the arrays so that they may be passed to "axes.bar3d".
    xData = xData.flatten()
    yData = yData.flatten()
    zData = zData.flatten()

    # Create figure
    fig = plt.figure()
    axes = fig.gca(projection='3d')
    axes.bar3d( xData, yData, np.zeros(len(zData)), .98, .98, zData, alpha=1.0, zsort='min')
      
    # Axes  
    axes.set_xlim3d(0, width)
    axes.set_ylim3d(0, height)
    
    # View
    if axesView[0] != 0 or axesView[1] != 0:
        axes.view_init(axesView[0], axesView[1])
    
    # Shift the tick labels up by minimum, set_zlim3d does not work
    if zRange[0] != 0 or zRange[1] != 0:
        if zTicks:
            axes.zaxis.set_major_formatter(ticker.FuncFormatter(major_formatter))
        else:
            axes.zaxis.set_major_formatter(ticker.NullFormatter())
   
    # Set the shifted range
    if zRange[0] != 0 or zRange[1] != 0:
        axes.set_zlim3d([0, zRange[1] - zRange[0]])

    plt.show()

# Plot a surface of data
def plotColorSurface(dataZ, colorsRGB, zRange = [0,0], stride = 1): 
    
    # Size of data edges. Edges is one more than matches
    width = dataZ.shape[1]
    height = dataZ.shape[0]
    
    # The edges of the surface. One more than the patches
    surfaceEdges = zeros((height + 1, width + 1), dtype='float')
       
    # Create  edge data form the height data
    surfaceEdges[0, 0] = dataZ[0, 0]
    surfaceEdges[0, width] = dataZ[0, width-1]
    surfaceEdges[height, 0] = dataZ[height-1, 0]
    surfaceEdges[height, width] = dataZ[height-1, width-1]
     
    for x in range(1, width):
        surfaceEdges[0, x] = (dataZ[0, x-1] + dataZ[0, x]) / 2.0
        surfaceEdges[height, x] = (dataZ[height-1, x-1] + dataZ[height-1, x]) / 2.0

    for y in range(1, height):
        surfaceEdges[y, 0] = (dataZ[y-1, 0] + dataZ[y, 0]) / 2.0
        surfaceEdges[y, width] = (dataZ[y-1, width-1] + dataZ[y, width-1]) / 2.0
       
    for x in range(1, width):
        for y in range(1, height):
            surfaceEdges[y, x] = (dataZ[y-1,x-1] + dataZ[y-1,x] + dataZ[y,x-1] + dataZ[y,x]) / 4.0

    # Create the x and y pixel indices arrays 
    x = np.linspace(0, width, width+1)
    y = np.linspace(0, height, height+1)
    xv, yv = np.meshgrid(x, y)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca(projection='3d')
    axes.plot_surface(xv, yv, surfaceEdges, rstride = stride, cstride = stride, linewidth=1, facecolors=colorsRGB, alpha=1.0, antialiased=False)
    
    axes.set_xlim3d(0, width)
    axes.set_ylim3d(0, height)
    if zRange[0] != 0 or zRange[1] != 0:
        axes.set_zlim3d([zRange[0], zRange[1]])

    plt.show()

# Plot a surface of data
def plotSurface(dataZ, zRange = [0,0], stride = 1): 
    
    # Size of data edges. Edges is one more than matches
    width = dataZ.shape[1]
    height = dataZ.shape[0]
    
    # The edges of the surface. One more than the patches
    surfaceEdges = zeros((height + 1, width + 1), dtype='float')
       
    # Create  edge data form the heigh data
    surfaceEdges[0, 0] = dataZ[0, 0]
    surfaceEdges[0, width] = dataZ[0, width-1]
    surfaceEdges[height, 0] = dataZ[height-1, 0]
    surfaceEdges[height, width] = dataZ[height-1, width-1]
     
    for x in range(1, width):
        surfaceEdges[0, x] = (dataZ[0, x-1] + dataZ[0, x]) / 2.0
        surfaceEdges[height, x] = (dataZ[height-1, x-1] + dataZ[height-1, x]) / 2.0

    for y in range(1, height):
        surfaceEdges[y, 0] = (dataZ[y-1, 0] + dataZ[y, 0]) / 2.0
        surfaceEdges[y, width] = (dataZ[y-1, width-1] + dataZ[y, width-1]) / 2.0
       
    for x in range(1, width):
        for y in range(1, height):
            surfaceEdges[y, x] = (dataZ[y-1,x-1] + dataZ[y-1,x] + dataZ[y,x-1] + dataZ[y,x]) / 4.0

    # Create the x and y pixel indices arrays 
    x = np.linspace(0, width, width+1)
    y = np.linspace(0, height, height+1)
    xv, yv = np.meshgrid(x, y)
    
    # Create figure
    fig = plt.figure()
    axes = fig.gca(projection='3d')
    axes.plot_surface(xv, yv, surfaceEdges, rstride = stride, cstride = stride, linewidth=1, alpha=1.0, antialiased=False)
    
    axes.set_xlim3d(0, width)
    axes.set_ylim3d(0, height)
    if zRange[0] != 0 or zRange[1] != 0:
        axes.set_zlim3d([zRange[0], zRange[1]])

    plt.show()
    
# Plot a surface of data
def plotWireframe(dataZ, zRange = [0, 0], stride = 1): 
    
    # Size of data edges. Edges is one more than matches
    width = dataZ.shape[1]
    height = dataZ.shape[0]
    
    # Create an X-Y mesh of the same dimension as the 2D data
    xData, yData = np.meshgrid( np.arange(width), np.arange(height))
    
    fig = plt.figure()    
    
    axes = fig.gca(projection='3d')
    axes.plot_wireframe(xData, yData, dataZ, rstride=stride, cstride=stride)
    
    axes.set_xlim3d(0, width)
    axes.set_ylim3d(0, height)
    if zRange[0] != 0 or zRange[1] != 0:
        axes.set_zlim3d([zRange[0], zRange[1]])

    plt.show()
    
# Plot a surface of data
def plotQuiver(magnitude, direction, scaleVectors = 0, sampleSpace = 1):
        
    width = magnitude.shape[1]
    height = magnitude.shape[0]                 
    
    # Create an X-Y mesh of the same dimension as the 2D data
    xPos, yPos = np.meshgrid( np.arange(width), np.arange(height))
     
    u = zeros((height, width), dtype='float')    
    v = zeros((height, width), dtype='float')    
    for y in range(0, height):
        for x in range(0, width):  
            # Image origin is top left and graph origin is bottom left, so invert and width units
            u[height - y - 1,x] = magnitude[y,x] * math.cos(direction[y,x]) /width
            v[height - y - 1,x] = magnitude[y,x] * math.sin(direction[y,x]) /width
    
    fig = plt.figure()
    axes = fig.gca()
    
    # For the upside down image
    initY = (height-1) % sampleSpace
     

    plt.quiver(xPos[initY::sampleSpace, ::sampleSpace], yPos[initY::sampleSpace, ::sampleSpace],   \
                u[initY::sampleSpace, ::sampleSpace], v[initY::sampleSpace, ::sampleSpace],        \
                pivot='mid', units='width', scale = scaleVectors)                                 
         

    axes.set_xlim(0, width-1)
    axes.set_ylim(0, height-1)    
    plt.show() 
