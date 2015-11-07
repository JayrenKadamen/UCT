__author__ = 'KDMJUN001'

import numpy as np
from calculations import leastSquares
from plotting import Display
import math as math
def main():
    #Getting provisional values (x0,y0,z0,r0)
    xlist = []
    ylist = []
    zlist = []

    file = open('data_assignment1.csv','r')

    for line in file:
        line = line.split(',')
        x_row = np.float(line[0])
        y_row = np.float(line[1])
        z_row = np.float(line[2])

        xlist+=[x_row]
        ylist+=[y_row]
        zlist+=[z_row]

    minx,miny,minz = min(xlist),min(ylist),min(zlist)
    maxx,maxy,maxz = max(xlist),max(ylist),max(zlist)

    x0 = (maxx-minx)/2 + minx
    y0 = (maxy-miny)/2 + miny
    z0 =  minz
    r0 = ((maxx-minx)/2 + (maxy-miny)/2+(maxz-minz))/3

    count = 1
    file.close()
    plotting = Display(xlist, ylist, zlist)

    iterate(x0,y0,z0,r0,count)
    plotting.plot()
def iterate(x0,y0,z0,r0,count):

    provs = leastSquares(x0,y0,z0,r0)
    x, V, A, M, B, B_Col  = provs.calc()

    prov = np.matrix([[x0],[y0],[z0],[r0]])

    soln = prov + x

    new_x0 = soln.item((0,0))
    new_y0 = soln.item((1,0))
    new_z0 = soln.item((2,0))
    new_r0 = soln.item((3,0))


    if abs(x.item((0,0))) and abs(x.item((1,0))) and abs(x.item((2,0))) and abs(x.item((3,0))) > .005:

        if count < 100:
            count += 1
            iterate(new_x0,new_y0,new_z0,new_r0,count)
    else:
        #Cofactor Matrix
        Qxx = ((A.T)*(M.I)*A).I
        S = (M.I)*A*Qxx
        Qkk = M.I - (M.I)*A*(S.T)

        #Reference Variance
        Var = ((V.T)*V)/(B_Col-3)
        #Covariance matrix of discrepancies w
        Eww = (Var.item((0,0)))*B*B.T
        #Covariance of Unknowns
        Exx = (Var.item((0,0)))*Qxx

        SigmaR = math.sqrt(Exx.item((3,3)))


        print Var, "Variance\n"
        print Exx, "Variance-Covariance of the Unknowns\n"
        print Eww, "Variance-Covariance of discrepancies\n"

        print count, "Iterations\n"
        print new_x0
        print new_y0
        print new_z0
        print new_r0
        provs.stats(SigmaR, B_Col)
main()