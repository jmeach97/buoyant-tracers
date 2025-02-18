# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:13:20 2022

@author: Jamie
"""
import numpy as np

for i in [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]:
    x,y=np.polynomial.legendre.leggauss(i)
    b=1
    a=0

    x=(b-a)*x/2+(a+b)/2
    y=y*(b-a)/2
    print(x)
    print(y)

    np.savetxt('GaussLegPoints'+str(i)+'.dat',x,delimiter=',')
    np.savetxt('GaussLegWeights'+str(i)+'.dat',y,delimiter=',')