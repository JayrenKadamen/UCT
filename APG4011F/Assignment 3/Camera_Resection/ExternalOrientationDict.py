'''
Created on 01 May 2015

@author: KDMJUN001
'''
from ExternalOrientation import Exterior


class ExternalOrientation(dict):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def keys(self,array):
        #x0,y0,z0,omega,phi,kappa,scale,name
        name = array[7]
        U0 = array[0]
        V0 = array[1]
        W0 = array[2]
        
        omega = array[3]
        phi = array[4]
        kappa = array[5]
        
        scale = array[6]
        
        self[name] = Exterior(U0,V0,W0,phi,omega,kappa,scale)