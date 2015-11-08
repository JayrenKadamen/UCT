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

names = ['Postgrad Lab','EGS','Snape','Average']
values = [0.997,1.027,1.015,1.013]

if __name__ == '__main__':

 
    trace1 = go.Bar(
        x=names,
        y=values,
marker=dict(
        color=['rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(222,45,38,0.8)'],
    ),
    )

    data = [trace1]
    layout = go.Layout(
                       title='Bar Graph of Scales per Room',
        font=dict(family='Raleway, sans-serif'),
        annotations=[
                dict(
                    x=xi,
                    y=yi,
                    text=str(yi),
                    xanchor='center',
                    yanchor='bottom',
                    showarrow=False,
                ) for xi, yi in zip(names, values)],
        xaxis=dict(
                   title='Scene'),
        yaxis=dict(
                   title='Scale',
                   range=[0.850,1.05],
                    titlefont=dict(
                                    size=18,
                                    color='rgb(107, 107, 107)'
                                    ),
                   showticklabels=True,
                   showgrid=False
                   ),

#         bargap=0.10,
#         barmode='group',
#         bargroupgap=0.0
    )
    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='scales')
    
    