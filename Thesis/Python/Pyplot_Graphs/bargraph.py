'''
Created on 25 Oct 2015

@author: jayre
'''
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api
import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np
from matplotlib import pyplot as plt
names = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']
ave = []
weights = []
labels = []

ave1 = []
weights1 = []
labels1 = []

ave2 = []
weights2 = []
labels2 = []

Bins = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1])

np.set_printoptions(suppress=False)

if __name__ == '__main__':

    folder1 = open("Deep_Space/H_Plane_data.csv")
#     folder2 = open("GTL/H_Plane_data.csv")
    folder2 = open("Snape/H_Plane_data.csv")
    folder3 = open("EGS/H_Plane_data.csv")
    
#     data = folder1
    
    for line in folder1:
        item = line.split(";")

        weights += [float(item[2])]
        
        ave += [float(item[10])]
        
        labels += [str(item[9])]
    
    ave_Array = np.array(ave)
    weights_Array = np.array(weights)
    
    hist, bin_edges = np.histogram(a=ave_Array, bins = Bins, weights = weights_Array) 
        
    for line in folder2:
        item = line.split(";")

        weights1 += [float(item[2])]
        
        ave1 += [float(item[10])]
        
        labels1 += [str(item[9])]
        
    ave_Array1 = np.array(ave1)
    weights_Array1 = np.array(weights1)
    
    hist1, bin_edges1 = np.histogram(a=ave_Array1, bins = Bins, weights = weights_Array1) 
    
    for line in folder3:
        item = line.split(";")

        weights2 += [float(item[2])]
        
        ave2 += [float(item[10])]
        
        labels2 += [str(item[9])]
    
    
    ave_Array2 = np.array(ave2)
    weights_Array2 = np.array(weights2)
    
    hist2, bin_edges2 = np.histogram(a=ave_Array2, bins = Bins, weights = weights_Array2)   
#-------------------------------------------------------------------------------------------------------------------

    
    trace1 = go.Bar(
        x=bin_edges,
        y=hist,
        name='Postgrad Lab',
        opacity=1
    )
    trace2 = go.Bar(
        x=bin_edges1,
        y=hist1,
        name='Snape',
        opacity=1
    )
    trace3 = go.Bar(
        x=bin_edges2,
        y=hist2,
        name='EGS',
        opacity=1
    )
    data = [trace1, trace2,trace3]
    layout = go.Layout(
                       title='Weighted Histogram of Horizontal Height Segments per Room',
        font=dict(family='Raleway, sans-serif'),
        xaxis=dict(
                   title='Height in Metres',
                   gridwidth=20,
                   autorange=True,
                   dtick=0.25,
                   tickangle=-45
                   ),
        yaxis=dict(
                   title='Point Count per Segment',
                    titlefont=dict(
                                    size=18,
                                    color='rgb(107, 107, 107)'
                                    ),
                   showticklabels=True
                   ),
        bargap=0.10,
        barmode='group',
        bargroupgap=0.0
    )
    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='grouped-bar')
    
    