'''
Created on 16 Sep 2015

@author: jayre
'''
import numpy as np
from matplotlib import pyplot as plt
ave = []

minH = []
minH_W = []

maxH = []
maxH_W = []


minV = []
minV_W = []

maxV = []
maxV_W = []

labels = []

Bins = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0])

np.set_printoptions(suppress=False)

if __name__ == '__main__':

    folder1 = open("Deep_Space/H_Plane_data.csv")
    folder2 = open("Snape/H_Plane_data.csv")
    folder3 = open("EGS/H_Plane_data.csv")
    
    data = folder1      
    
    for line in data:
        item = line.split(",")

        minH += [float(item[6])]
        minH_W += [float(item[2])]
        
        maxH += [float(item[8])]
        maxH_W += [float(item[2])]
        
        ave += [float(item[10])]
        
        labels += [str(item[9])]
    
    maxH += [float(item[8])]
    maxH_W += [float(item[2])]
    
    minH_Array = np.array(minH)
    minH_W_Array = np.array(minH_W)
    
    maxH_Array = np.array(maxH)
    maxH_W_Array = np.array(maxH_W)
    
    ave_Array = np.array(ave)
    
    hist, bin_edges = np.histogram(a=minH_Array, bins = Bins, weights = minH_W_Array)
    hist2,bin_edges2 = np.histogram(a=maxH_Array, bins = Bins, weights = maxH_W_Array)
    hist3,bin_edges3 = np.histogram(a=ave_Array, bins = Bins, weights = maxH_W_Array)
    
#     Plotting Min and Max Seperately
    f = plt.figure(figsize = plt.figaspect(0.45))
    ax = f.add_subplot(1,2,1)
    plt.title("Minimum Heights")
    ax.bar(bin_edges[:-1],hist,width=0.1)
    ax.set_xlim(0,max(bin_edges))
#     ax.set_ylim([min(hist),max(hist)+10000])
     
    ax = f.add_subplot(1,2,2)
    plt.title("Maximum Heights")
    ax.bar(bin_edges2[:-1],hist2,width=0.1)
    ax.set_xlim(0,max(bin_edges2))
#     ax.set_ylim([min(hist),max(hist)+10000])
    
#     f.show()

 
    combined_H = np.append(arr=minH_Array,values=maxH_Array)
    combined_W = np.append(arr=minH_W_Array, values=maxH_W_Array)
    
    histC, bin_edgesC = np.histogram(a=combined_H, bins = Bins, weights = combined_W)
    
    
    #Plotting Min and Max on same axes
    g = plt.figure(2)
    plt.title("Combined Histogram of Heights")
    plt.bar(bin_edgesC[:-1],histC,width=0.1)
    plt.xlim(0,max(bin_edges))
    
    
#     g.show()

    g = plt.figure(3)
    plt.title("Histogram of Average Heights")
    plt.bar(bin_edges3[:-1],hist3,width=0.1)
    plt.xlim(0,max(bin_edges3))
    
    plt.show()
    


    
            