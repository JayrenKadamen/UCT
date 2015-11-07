'''
Created on 28 Apr 2015

@author: KDMJUN001
'''


class Exterior(object):
    '''
    classdocs
    '''
    

    def __init__(self, U0,V0,W0,phi,omega,kappa,scale):
        '''
        Constructor
        '''
        self.U0 = U0
        self.V0 = V0
        self.W0 = W0
        
        self.omega = omega
        self.phi = phi
        self.kappa = kappa
    
        self.scale = scale
        
