'''
Created on 01 May 2015

@author: KDMJUN001
'''
from InternalGeometry import Interior


class InternalOrientation(dict):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def keys(self,array):
        #name,c,py,px,pxl_size
        name = array[0]
        c = array[1]
        
        #Image size (physical size)
        py = array[2]
        px = array[3]

        
        self[name] = Interior(c,py,px)