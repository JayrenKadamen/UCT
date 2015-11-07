'''
Created on 09 Apr 2015

@author: KDMJUN001
'''
import numpy as np


class LeastSquare(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

        
    def normal(self,x,y,z):

        
        Alist = []
        Llist = []
        for i in range(0,len(x)):
            
            Alist += [[x[i],y[i], 1.0]]
            
            Llist += [[-z[i]]]
        
        A = np.matrix(Alist)
        
        L = np.matrix(Llist)
        
        #Return non singular matrix results
        if np.linalg.det(A.T*A) != 0: 
            pp = ((A.T * A).I * A.T * L)

            pp = pp.tolist()

            return [pp[0][0], pp[1][0], 1.0], pp[2][0]
        

        # If the matrix is singular then return a zero vector
        else: 

            return [0,0,0],0