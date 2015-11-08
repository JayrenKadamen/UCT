'''
Created on 16 Sep 2015

@author: jayre
'''
import numpy as np
from matplotlib import pyplot as plt
ave = []

weights = []

labels = []

Bins = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0])

np.set_printoptions(suppress=False)

if __name__ == '__main__':

    folder1 = open("Deep_Space/H_Plane_data.csv")
    folder2 = open("Snape/H_Plane_data.csv")
    folder3 = open("EGS/H_Plane_data.csv")
    
    data = folder3
    
    for line in data:
        item = line.split(";")

        weights += [float(item[2])]
        
        ave += [float(item[10])]
        
        labels += [str(item[9])]
    
    
    ave_Array = np.array(ave)
    weights_Array = np.array(weights)
    
    hist, bin_edges = np.histogram(a=ave_Array, bins = Bins, weights = weights_Array)    
    print(hist)
    print(bin_edges)
    #Plotting Min and Max on same axes
    g = plt.figure(1)
    plt.title("Histogram of Average Heights")
    plt.bar(bin_edges[:-1],hist,width=0.1)
    plt.xlim(0,max(bin_edges))
    
    
    plt.show()
    


    
            