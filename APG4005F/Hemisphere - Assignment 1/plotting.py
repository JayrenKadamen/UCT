__author__ = 'KDMJUN001'
'''
Created on 04 Mar 2015

@author: kdmjun001
'''

from matplotlib import pylab as plt
from matplotlib import cm as CM
from mpl_toolkits.mplot3d import axes3d


class Display(object):
    '''
    classdocs
    '''


    def __init__(self, xlist,ylist,zlist):

        self.X = xlist
        self.Y = ylist
        self.Z = zlist

    def plot(self):



        fig = plt.figure()

        ax = fig.add_subplot(1,1,1, projection='3d')
        plt.title('Scatter Plot of Measured Points')


        ax.scatter(self.X,self.Y,self.Z,s=50,c=self.Z,cmap=CM.jet)

        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')

        plt.show()


