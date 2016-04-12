# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\umesh123\.spyder2\.temp.py
"""
from scipy import linspace
import scipy
from pylab import plot,show,pause

'''A + B -- C'''      #the reaction considered
L=100.0               #Length of pfr
u=0.02                #flowrate in pfr
Da=2.0                #dispersion coefficient
k=0.2                 #rate constant
Ca0=1.0               #initial concenteration

t=0
dt=0.1
dz=L/100
Db=0.2

'''For nonideal PFR'''
C=scipy.ones(10)*0    #concenteartion of A
C[0]=1
C[1]=1

Cb=scipy.ones(10)*0   #concenteration of B
Cb[0]=1
Cb[1]=1

'''For ideal pfr'''   #just to compare with non ideal pfr 
Ci=scipy.ones(10)*0   #ideal concenteration of A
Ci[0]=1
Ci[1]=1

Cbi=scipy.ones(10)*0  #ideal concenteration of B
Cbi[0]=1
Cbi[1]=1

l=linspace(0,L,10)
while t<10.0:
    for i in range(1,9):
        Ci[i+1]=Ci[i]-u*dt*(Ci[i]-Ci[i-1])/dz-k*dt*Ci[i]*Cbi[i]
        Cbi[i+1]=Cbi[i]-dt*u*(Cbi[i]-Cbi[i-1])/dz-k*dt*Ci[i]*Cbi[i]
        C[i+1]=C[i]+dt*Da*(C[i+1]-2*C[i]+C[i-1])/dz/dz-dt*u*(C[i]-C[i-1])/dz-k*dt*C[i]*Cb[i]
        Cb[i+1]=Cb[i]+dt*Db*(Cb[i+1]-2*Cb[i]+Cb[i-1])/dz/dz-dt*u*(Cb[i]-Cb[i-1])/dz-k*dt*C[i]*Cb[i]
    
    print ('The concenteration(ideal) of A is', Ci[i+1])
    print ('The concenteration(ideal) of B is', Cbi[i+1])
    print ('The concenteration(non ideal) of A is', C[i+1])
    print('The concenteration(non ideal) of B is', Cb[i+1])
    print ""
    plot(l,C)
    plot(l,Ci)
    show()
    pause(0.1)        #showing graphs for all the iterations
    
