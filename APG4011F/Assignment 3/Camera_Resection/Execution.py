'''
Created on 29 Apr 2015

@author: KDMJUN001
'''
#Importing created classes
import InternalGeometry as IG
from ExternalOrientationDict import ExternalOrientation
from InternalOrientationDict import InternalOrientation

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

import pylab


#Python Libraries
import numpy as np
from random import randint, random, uniform
import sympy as sy

sym_x,sym_y,sym_x0,sym_y0,sym_X0,sym_Y0,sym_Z0,sym_omega,sym_phi,sym_kappa,sym_c = sy.symbols('x,y,x0,y0,X0,Y0,Z0,omega,phi,kappa,c')

class Execute(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
    
    def Setup(self):
        
        self.imageGrid()
            
            
    def imageGrid(self):
        #Creation of n random points
        image1 = []
        global n
        n = 30
        for i in range(0,n):
            image1 += [[uniform(-0.3,0.3),uniform(-0.3,0.3)]]
        
        #Creation of scale for each point
        for i in range(0,n):
            image1[i].append((randint(90,100)))
        
        #Converting list to numpy array to allow for easy indexing
        global image1Array

        image1Array =  (np.asarray(image1))
        
#         print("Randomly Generated Image Points with Scale")
#         for i in range(0,n):
#             print (image1[i])
#         print()
        
        self.cameraInterior()
        
    def cameraInterior(self):
        global cameraInternal
        cameraInternal = InternalOrientation()
        
        #All values sent are in mm 
        #name, c,xi,yi
        IP = ['Canon',0.3,0,0]
        cameraInternal.keys(IP)
        
        self.cameraOrientation()
    
    def cameraOrientation(self):
        #Setting camera parameters for image 1
        global cameraExternal
        cameraExternal = ExternalOrientation()

        #Coordinates for Kings Battery = 50750.110, 58299.800
        #U0,V0,W0,phi,omega,kappa,scale,camera
        cameraOne = [0, 0, 10, np.radians(0),np.radians(0),np.radians(0),(1/100),'cameraOne']
        cameraExternal.keys(cameraOne)

        cameraTwo = [10, 0, 10,np.radians(0.5),np.radians(0.5),np.radians(0.05),(1/100),'cameraTwo']
        cameraExternal.keys(cameraTwo)
        
        self.ObjectPoints()
        
    
    def ObjectPoints(self):
        omega,phi,kappa= sy.symbols('omega,phi,kappa')
        
        omega = cameraExternal['cameraOne'].omega
        phi = cameraExternal['cameraOne'].phi
        kappa = cameraExternal['cameraOne'].kappa
 
        Rphi = np.matrix([[sy.cos(phi),0,-sy.sin(phi)],
                          [0,          1,           0],
                          [sy.sin(phi),0,sy.cos(phi)]])
        
        Romega = np.matrix([[1,            0,            0],
                            [0,sy.cos(omega), sy.sin(omega)],
                            [0,-sy.sin(omega),sy.cos(omega)]])
        
        Rkappa = np.matrix([[sy.cos(kappa),sy.sin(kappa), 0],
                            [-sy.sin(kappa),sy.cos(kappa),0],
                            [0,           0,              1]])
        global R
#         R = Romega*Rkappa*Rphi
        R = Rphi*Rkappa*Romega
        projectionCentre = [[cameraExternal['cameraOne'].U0],[cameraExternal['cameraOne'].V0],[cameraExternal['cameraOne'].W0]]
        
        global objectPoints1
        objectPoints1 = []
        
        for i in range(0,n):
            
            internal = [[image1Array[i,0]], #x
                        [image1Array[i,1]], #y
                        [-cameraInternal['Canon'].c]]   #focal length         
            
            point = image1Array[i,2]*R*internal+projectionCentre

            objectPoints1 += [point]
        
        objectPoints1 = np.asarray(objectPoints1)
        
#         for i in range (0,n):
#             print(objectPoints1[i][0],objectPoints1[i][1],objectPoints1[i][2],"\n")
        global objectPointsError
        objectPointsError = []
        
        #Adds small error to points
        for i in range(0,n):
            
            error = [objectPoints1[i,0]+(random()/10),objectPoints1[i,1]+(random()/10),objectPoints1[i,2]+(random()/10)]
            objectPointsError += [error]
            
        objectPointsError = np.asarray(objectPointsError)
        
        self.ImagePoints2()
        
    def ImagePoints2(self):
        
        image2 = []
        
        for i in range (0,n):

            x = (-cameraInternal['Canon'].c)*((R[0,0]*(objectPoints1[i,0]-cameraExternal['cameraTwo'].U0)
                 +R[0,1]*(objectPoints1[i,1]-cameraExternal['cameraTwo'].V0)
                 +R[0,2]*(objectPoints1[i,2]-cameraExternal['cameraTwo'].W0))/(R[2,0]*(objectPoints1[i,0]-cameraExternal['cameraTwo'].U0)+
                  R[2,1]*(objectPoints1[i,1]-cameraExternal['cameraTwo'].V0)+
                  R[2,2]*(objectPoints1[i,2]-cameraExternal['cameraTwo'].W0)))
            
            y = (-cameraInternal['Canon'].c)*((R[1,0]*(objectPoints1[i,0]-cameraExternal['cameraTwo'].U0)
                 +R[1,1]*(objectPoints1[i,1]-cameraExternal['cameraTwo'].V0)
                 +R[1,2]*(objectPoints1[i,2]-cameraExternal['cameraTwo'].W0))/(R[2,0]*(objectPoints1[i,0]-cameraExternal['cameraTwo'].U0)+
                  R[2,1]*(objectPoints1[i,1]-cameraExternal['cameraTwo'].V0)+
                  R[2,2]*(objectPoints1[i,2]-cameraExternal['cameraTwo'].W0)))
            
            image2 += [[x,y]]
            
        global image2Array
        image2Array = np.asarray(image2)
           
        global image2Error
        image2Error = []
        
        #Adds small error to points
        for i in range(0,n):
            
            error = [image2Array[i,0],image2Array[i,1]]
            image2Error += [error]
            
        image2Error = np.asarray(image2Error)        
        
        self.lies()
        self.output()
        
    def ObjectPoints2(self):
        projectionCentre = [[cameraExternal['cameraTwo'].U0],[cameraExternal['cameraTwo'].V0],[cameraExternal['cameraTwo'].W0]]
        
        global objectPoints2
        objectPoints2 = []
        
        for i in range(0,n):
            
            internal = [[image2Array[i,0]], #x
                        [image2Array[i,1]], #y
                        [-cameraInternal['Canon'].c]]   #focal length         
            
            point = image2Array[i,2]*R*internal+projectionCentre

            objectPoints2 += [point]
        
        objectPoints2 = np.asarray(objectPoints2)
        
        self.output()
        
    def output(self):
        fig = plt.figure(figsize=(80,60))
        ax = fig.add_subplot(1,1,1,projection ='3d')
        plt.title('Projection Centre (0,0,10), Principal Distance = 0.3m, Scale = 1:100')
        
        ax.scatter(0, 0, 10, c='blue')
        ax.scatter(10, 0, 10, c='red')

        
        for i in range(0,n):
#             #Image 1 points
#             ax.scatter(image1Array[i][0]+10,image1Array[i][1]+10, c='blue')
#             #Image 2 points
#             ax.scatter(image2Error[i][0]+10,image2Error[i][1], c='red')
            #Object points 1
            ax.scatter(objectPoints1[i][0],objectPoints1[i][1],objectPoints1[i][2], c='blue')
            #Object points 2
            ax.scatter(objectPointsError[i][0],objectPointsError[i][1],objectPointsError[i][2], c='red')            
            #Projection Rays - Image 1
            ax.plot([0,objectPoints1[i][0]],[0,objectPoints1[i][1]],[10,objectPoints1[i][2]], c='blue')
#             #Projection rays - Image 2
            ax.plot([10,objectPointsError[i][0]],[0,objectPointsError[i][1]],[10,objectPointsError[i][2]], c='red')
                    

#         plt.show()
        self.imagepointsplot()

    def imagepointsplot(self):
#        Code for 3D plot of image points
#         fig2 = plt.figure(figsize=(80,60))
#         bx = fig2.add_subplot(1,1,1,projection='3d')
        
#         for i in range(0,n):
#             bx.scatter(image1Array[i][0],image1Array[i][1], c='blue', s = 50)
#             bx.scatter(image2Error[i][0],image2Error[i][1], c='red', s = 50)

        for i in range(0,n):
            plt.scatter(image1Array[i][0],image1Array[i][1], c='green')
            plt.scatter(image2Error[i][0],image2Error[i][1], c='red')
            
#         plt.show()
        
    def lies(self):
        sum1 = []
        for i in range(0,n):
            one = objectPoints1[i][0]-objectPointsError[i][0]
            two = objectPoints1[i][1]-objectPointsError[i][1]
            three = objectPoints1[i][2]-objectPointsError[i][2]
            
            print(image1Array[i][0])
            print(image1Array[i][1])
            print(image2Error[i][0])
            print(image2Error[i][1],"\n")
            print (objectPoints1[i][0])
            print (objectPoints1[i][1])
            print (objectPoints1[i][2], "\n\n")
            sum1 += [one+two+three]
            
        mean_error = (sum(sum1)/(len(sum1)))
        sum2 = np.asarray(sum1)
        sum_sss = []
        
        for i in range(0,len(sum1)):
            sum_sss += [(sum2[i]-mean_error)]
            
        print(mean_error)
        print(1/(len(sum1))*(sum(sum_sss))**2)
        print()
        
        count = 0
        
        def errors(count,mean_error):
            x = random()
            if x < abs(mean_error):
                print (x)
                count += 1
                if count < 6:
                    errors(count,mean_error)
                elif count == 6:
                    pass
            else:
                errors(count,mean_error)
        errors(count,mean_error)
