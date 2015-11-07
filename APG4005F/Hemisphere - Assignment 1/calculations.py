__author__ = 'KDMJUN001'

import numpy as np
from scipy.stats import chi2


class leastSquares(object):
    '''
    classdocs
    '''

    def __init__(self,x0,y0,z0,r0):
        '''
        Constructor
        '''
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r0 = r0

    def calc(self):
        count = 0
        file = open('data_assignment1.csv','r')
        listA = []
        listB = []
        listw = []

        for line in file:
            line = line.split(',')
            x = np.float(line[0])
            y = np.float(line[1])
            z = np.float(line[2])

            dx0 = -2*(x-self.x0)
            dy0 = -2*(y-self.y0)
            dz0 = -2*(z-self.z0)
            dr0 = -2*self.r0

            dx = 2*(x-self.x0)
            dy = 2*(y-self.y0)
            dz = 2*(z-self.z0)

            fn = (x-self.x0)**2+(y-self.y0)**2+ (z-self.z0)**2-(self.r0)**2

            listA += [[dx0,dy0,dz0,dr0]]
            listB += [[dx,dy,dz]]
            listw += [[fn]]

        file.close()
        B_Col = len(listA)*3

        BArray = []

        index = 0
        for row in listB:
            rowArray = [0]*B_Col
            rowArray[index+0] = row[0]
            rowArray[index+1] = row[1]
            rowArray[index+2] = row[2]
            BArray+= [rowArray]
            index +=3

        A = np.matrix(listA)
        B = np.matrix(BArray)
        W = np.matrix(listw)

        M = B*(B.T)

        X = -(A.T*(M.I)*A).I*A.T*M.I*W

        #Vector of correlates (k)
        k = -(M.I)*(A*X + W)
        #Vector of residuals
        V = (B.T)*k

        return X, V, A, M, B, B_Col
    
    def stats(self,SigmaR,B_Col):
        #Chi Squared test
        print('Our Null Hypothesis states that the variance of our population = sample variance')

        #Variance of our radius from measured points
        observed = SigmaR

        #Expected variance (3mm per x,y,z observation)
        expected = .003**2+.003**2+.003**2

		#Calculation of degrees of freedom
        dof = B_Col-1

        #Calculation of test statistics
        teststatx = B_Col*((observed - expected)**2/expected)
        teststatx1 = dof*(observed/expected)

        #User is prompted to input desired significance level
        significance = np.float(input('Please specify the significance level: '))

        print(teststatx),
        print(teststatx1)

        #Using built in scipy.stats.chi2 function instead of looking up values on a table
        mean, var, skew, kurt = chi2.stats(dof, moments='mvsk')
        Chi = chi2.ppf((1-significance),dof)

		#If our sampled variance is greater than the population variance at the chosen significance level then we reject the null hypothesis at that significance level
        if teststatx > Chi:
            print 'We reject the null hypothesis at the ',significance,'significance level'

        else:
            print 'We fail to reject the null hypothesis at the ',significance,'significance level'

        print(teststatx, dof)

