'''
Created on Apr 9, 2015

@author: kdmjun001
'''
import numpy as np
from scipy.spatial import KDTree
from LeastSquares import LeastSquare
import datetime
#For diagnostic computing time only
import time


def main():

    start_time_1 = time.time()
    
    output = open('Results.csv','w')
    myfile = open('W55A - Cloud.csv','r')
    
    coordinates = []
    
    for line in myfile:
        element = line.split(',')
        x = float(element[0])
        y = float(element[1])
        z = float(element[2])

        #Creates 2D list of points from specified file
        coordinates += [[x,y,z]]
    
    import_time = (time.time()-start_time_1)
    
    #Allows user to track the progress of the program execution
    percent = 0
    
    #Creates an index which can be used to look up the nearest neighbour of any point in given array
    tree = KDTree(coordinates)
    
    start_time_2 = time.time()
    #Here we are taking each point in the coordinates array...
    for point in coordinates:
        #...and finding the 30 nearest neighbours to that point
        k = 30
        index = tree.query(point,k)[1]
        
        X,Y,Z = [],[],[]

        for i in list(index):

            #X,Y,Z are lists containing the index of the 10 nearest neighbours to the point currently selected in the loop
            X += [coordinates[i][0]]
            Y += [coordinates[i][1]]
            Z += [coordinates[i][2]]


        #Calculating the plane and the normal to the plane, please see LeastSquares Class for more information
        leastsquares = LeastSquare()
        newnormal = leastsquares.normal(X,Y,Z)[0]

        # If the LeastSquares function returns anything but a zero vector, calculate the normal to that vector
        if newnormal != ([0,0,0],0):
            
            normal = np.array(newnormal)

            #Returns one of an infinite number of vector norms
            vnorm = np.linalg.norm(normal)
            #Calculate unit vector
            normalunit = normal/vnorm

            vert = [0,0,1]
            
            #angle between unit_normal and vertical
            vert_ang = np.arccos( (normalunit[0] * vert[0]  + normalunit[1] * vert[1] + normalunit[2] * vert[2]) )
              
            #Write points,normal,and angle to file
            output.writelines(str(point[0])+','+str(point[1])+','+str(point[2])+','+str(normal[0]) +','+str(normal[1]) +','+str(normal[2]) +',' +str(np.rad2deg(vert_ang))+'\n')    
            
            #Display program progress as a percentage of processed points
            percent += (1.0/len(coordinates))*100.0
            percentRound = round(percent,3)
            
#             if percentRound.is_integer():
#                 print ("Progress: "+str(percentRound)+' %')
#             else:
#                 pass

            print ("Progress: "+str(percentRound)+' %')
        else:
            pass
    computation_time = (start_time_2-start_time_1)
    output.close() 
    
    print ('Complete')
    print (import_time," seconds for reading file of ", len(coordinates), "points")
    print (computation_time," seconds for performing nearest neighbour calculation of", k, "neighbours per point")
    
    timetaken = str(datetime.timedelta(seconds=((time.time()-start_time_1))))
    timestr = timetaken.split(':')
    hours = timestr[0]
    mins = timestr[1]
    seconds = round(float(timestr[2]),2)
    
    print (hours," Hours", mins," minutes",seconds," seconds")


main()




    




