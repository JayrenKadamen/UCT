'''
Created on 16 Sep 2015

@author: jayre
'''
import numpy as np
from matplotlib import pyplot as plt

Weight = []
Height = []
NormHeight = []

labels = []

Bins = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0])
NormBin = np.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])


np.set_printoptions(suppress=False)

if __name__ == '__main__':

    folder1 = open("Deep_Space/DeepSpace.csv")
    folder2 = open("EGS/EGS.csv")
    folder3 = open("Snape/Snape.csv")

    
    data = folder1
    data2 = folder2
    
    for line in data:
        item = line.split(";")
        
        Weight += [float(item[2])]
        Height += [float(item[10])]
        NormHeight += [float(item[12])]


        
        labels += [str(item[11])]
    
    
    Weight_Array = np.array(Weight)
    Height_Array = np.array(Height)    
    NormHeight_Array = np.array(NormHeight)

    hist, bin_edges = np.histogram(a=Height_Array, bins = Bins, weights = Weight_Array)
    histNorm,bin_edges_Norm = np.histogram(a=NormHeight, bins = NormBin, weights = Weight_Array)


    #Plotting Min and Max on same axes
    f = plt.figure(figsize = plt.figaspect(0.45))
    ax = f.add_subplot(1,2,1)
#     f = plt.figure(1)
    plt.title("Histogram of Average Segment Heights")
    plt.bar(bin_edges[:-1],hist,width=0.1)
    plt.xlim(0,max(bin_edges))
#     plt.ylim(0,max(hist)+20000)
    
    ax = f.add_subplot(1,2,2)
#     g = plt.figure(2)
    plt.title("Normalised Histogram of Average Segment Heights")
    plt.bar(bin_edges_Norm[:-1],histNorm,width=0.05)
    plt.xlim(0,max(bin_edges_Norm))
#     plt.ylim(0,max(hist)+20000)

    for line in data2:
        item = line.split(";")
        
        Weight += [float(item[2])]
        Height += [float(item[10])]
        NormHeight += [float(item[12])]


        
        labels += [str(item[11])]
    
    
    Weight_Array = np.array(Weight)
    Height_Array = np.array(Height)    
    NormHeight_Array = np.array(NormHeight)

    hist, bin_edges = np.histogram(a=Height_Array, bins = Bins, weights = Weight_Array)
    histNorm,bin_edges_Norm = np.histogram(a=NormHeight, bins = NormBin, weights = Weight_Array)


    #Plotting Min and Max on same axes
    g = plt.figure(figsize = plt.figaspect(0.45))
    ax = g.add_subplot(1,2,1)
#     f = plt.figure(1)
    plt.title("Histogram of Average Segment Heights")
    plt.bar(bin_edges[:-1],hist,width=0.1)
    plt.xlim(0,max(bin_edges))
#     plt.ylim(0,max(hist)+20000)
    
    ax = g.add_subplot(1,2,2)
#     g = plt.figure(2)
    plt.title("Normalised Histogram of Average Segment Heights")
    plt.bar(bin_edges_Norm[:-1],histNorm,width=0.05)
    plt.xlim(0,max(bin_edges_Norm))
#     plt.ylim(0,max(hist)+20000)
    
    plt.show()
    


    
            