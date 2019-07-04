'''
Feature Extraction and Image Processing 
Mark S. Nixon & Alberto S. Aguado
http://www.southampton.ac.uk/~msn/book/
GeometricUtilities: Helper module to apply geometric transformations to images
'''

# Images
from ImageUtilities import createImageRGB

# Math
from math import sqrt

# Iteration
from timeit import itertools

# Image data
import numpy as np


# Solves the equation A*x=b
def solveSystem(A, b):
    x = np.linalg.solve(A, b)
    return x

# Transform an image
def imageTransform(image, maskImage, T):    
    height, width = len(image), len(image[0])
    centreX, centreY = width/2, height/2
    
    sImage = createImageRGB(width, height)
    for y, x in itertools.product(range(0, height-1), range(0, width-1)):
        # Alpha and colour   
        alpha = maskImage[y,x]/256.0 
        if alpha == 0: 
            continue
        rgb = (image[y,x]/4.0 + image[y+1,x+1]/4.0 + image[y+1,x]/4.0 +  \
               image[y,x+1]/4.0) * alpha
        
        # Transform
        cx, cy = x - centreX, y - centreY
        p0z = T[2][0] * cx + T[2][1] * cy + T[2][2] 
        p1z = T[2][0] * (cx+1) + T[2][1] * cy + T[2][2] 
        p2z = T[2][0] * (cx+1) + T[2][1] * (cy+1) + T[2][2] 
        
        if p0z != 0 and p1z != 0 and p2z !=0:
    
            p0x = int((T[0][0] * cx + T[0][1] * cy + T[0][2]) / p0z + centreX)
            p0y = int((T[1][0] * cx + T[1][1] * cy + T[1][2]) / p0z + centreY) 
        
            p1x = int((T[0][0] * (cx+1) + T[0][1] * cy + T[0][2]) / p1z + centreX)
            p1y = int((T[1][0] * (cx+1) + T[1][1] * cy + T[1][2]) / p1z + centreY) 
        
            p2x = int((T[0][0] * (cx+1) + T[0][1] * (cy+1) + T[0][2]) / p2z + centreX)
            p2y = int((T[1][0] * (cx+1) + T[1][1] * (cy+1) + T[1][2]) / p2z + centreY) 
            
            # Fill output image
            v1, v2 = [p1x - p0x, p1y - p0y], [p2x - p0x, p2y - p0y]
                    
            lv1 = max(.001,sqrt(v1[0]*v1[0] + v1[1]*v1[1]))
            lv2 = max(.001,sqrt(v2[0]*v2[0] + v2[1]*v2[1]))
            v1N = [v1[0]/lv1, v1[1]/lv1]
            v2N = [v2[0]/lv2, v2[1]/lv2]
    
            for dV1, dV2 in itertools.product(range(0, int(lv1)+1), range(0, int(lv2)+1)):
                a,b = int(p0x + dV1 * v1N[0] + dV2 * v2N[0]), int(p0y + dV1 * v1N[1] + dV2 * v2N[1])
                if a>0 and a < width and b > 0 and b < height:
                    sImage[b,a] = rgb
    return sImage   


# Get corresponding points in the 3d points plane origin,v1,v2
def projectionPoints(origin, v1, v2, npts, p, centreX, centreY):
    xy = [ ]
    for a in range(0, npts):
        rowxy = [ ]
        for b in range(0, npts):
        
            v1D = [a*v1[0]/float(npts-1), a*v1[1]/float(npts-1), a*v1[2]/float(npts-1)]
            v2D = [b*v2[0]/float(npts-1), b*v2[1]/float(npts-1), b*v2[2]/float(npts-1)]
            s = [origin[0]+v1D[0]+v2D[0], origin[1]+v1D[1]+v2D[1], origin[2]+v1D[2]+v2D[2]]
        
            sx = p[0]*s[0] + p[1]*s[1] + p[2]*s[2]  + p[3]
            sy = p[4]*s[0] + p[5]*s[1] + p[6]*s[2]  + p[7]
            sz = p[8]*s[0] + p[9]*s[1] + p[10]*s[2] + p[11] 
 
            rowxy.append([int((sx/sz) + centreX), int((sy/sz) + centreY)]) 
   
        xy.append(rowxy)
    
    return xy

# Get all the projection points of a unit cube
def projectionCubePoints(npts, p, centreX, centreY):
    xy1 = projectionPoints([0,0,1],[1,0,0], [0,1,0], npts, p, centreX, centreY)
    xy2 = projectionPoints([0,1,0],[1,0,0], [0,0,1], npts, p, centreX, centreY) 
    xy3 = projectionPoints([1,0,0],[0,1,0], [0,0,1], npts, p, centreX, centreY)
    
    return [xy1,xy2,xy3]   
        
def fillImage(colour, xy, image):
    npts = len(xy)
    
    height = len(image)
    width = len(image[0]) 
    
    for a,b in itertools.product(range(0, npts-1), range(0, npts-1)):
        
        c0,c1,c2 = xy[a][b], xy[a+1][b], xy[a][b+1]   
        
        # Fill output image
        v1 = [c1[0]-c0[0], c1[1]-c0[1]]
        v2 = [c2[0]-c0[0], c2[1]-c0[1]]
        
        lv1 = max(.001,sqrt(v1[0]*v1[0] + v1[1]*v1[1]))
        lv2 = max(.001,sqrt(v2[0]*v2[0] + v2[1]*v2[1]))
        v1N = [v1[0]/lv1, v1[1]/lv1]
        v2N = [v2[0]/lv2, v2[1]/lv2]
        
        for dV1, dV2 in itertools.product(range(0, 4*(1+int(lv1))), range(0, 4*(1+int(lv2)))):
            x = int(c0[0] + v2N[0] * dV2*.25 + v1N[0] * dV1*.25)
            y = int(c0[1] + v2N[1] * dV2*.25 + v1N[1] * dV1*.25)
            if x>0 and x < width and y > 0 and y < height:
                image[y,x] = colour
                
         
def fillImageColours(colours, xy, image):
    nfaces = len(xy)
    
    height = len(image)
    width = len(image[0]) 
    
    for faceNum in range(0, nfaces):
        face = xy[faceNum]
        npts = len(face)
        colour = colours[faceNum]
    
        for a,b in itertools.product(range(0, npts-1), range(0, npts-1)):
            
            c0,c1,c2 = face[a][b], face[a+1][b], face[a][b+1]  

            c = colour[a,b];
            
            # Fill output image
            v1 = [c1[0]-c0[0], c1[1]-c0[1]]
            v2 = [c2[0]-c0[0], c2[1]-c0[1]]
            
            lv1 = max(.001,sqrt(v1[0]*v1[0] + v1[1]*v1[1]))
            lv2 = max(.001,sqrt(v2[0]*v2[0] + v2[1]*v2[1]))
            v1N = [v1[0]/lv1, v1[1]/lv1]
            v2N = [v2[0]/lv2, v2[1]/lv2]
            
            for dV1, dV2 in itertools.product(range(0, 4*(1+int(lv1))), range(0, 4*(1+int(lv2)))):
                x = int(c0[0] + v2N[0] * dV2*.25 + v1N[0] * dV1*.25)
                y = int(c0[1] + v2N[1] * dV2*.25 + v1N[1] * dV1*.25)
                if x>0 and x < width and y > 0 and y < height:
                    image[y,x] = c
        
def computeProjection(pts,q):
    # Fill matrix
    M = [ ]    
    for row in range(0,6):
        r1 = [ q[row][0],q[row][1], q[row][2],1,0,0,0,0,-pts[row][0]*q[row][0],  \
              -pts[row][0]*q[row][1],-pts[row][0]*q[row][2],-pts[row][0] ]
        r2 = [ 0,0,0,0,q[row][0],q[row][1], q[row][2],1,-pts[row][1]*q[row][0],  \
              -pts[row][1]*q[row][1],-pts[row][1]*q[row][2],-pts[row][1] ]
        M.append(r1)
        M.append(r2)
    print(M)

    # Solves the equation A*x=b
    r = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]
    p = solveSystem(M, r)
   
    return p

def getPointColours(xy, mask, image):
    nfaces = len(xy)
    
    height = len(image)
    width = len(image[0]) 
    
    colourImages = [ ]
    for faceNum in range(0, nfaces):
    
        face = xy[faceNum]
        npts = len(face)
        colours = createImageRGB(npts, npts)

        for a,b in itertools.product(range(0, npts-1), range(0, npts-1)):
            y,x = face[a][b][1], face[a][b][0]
            if y>0 and y<height and x>0 and x<width:
                alpha = mask[y,x] / 256.0 
                if alpha > 0.0:
                    c = image[y,x]
                    colours[a,b] = [alpha*c[0],alpha*c[1],alpha*c[2]]
                
        colourImages.append(colours)

    return colourImages
    
