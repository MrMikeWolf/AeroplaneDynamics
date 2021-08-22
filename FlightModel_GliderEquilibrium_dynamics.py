# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 21:04:49 2021

@author: H295476
"""
import numpy as np
from numpy import sin,cos
from scipy.integrate import odeint
import matplotlib.pyplot as plot
import numpy as np

g = 32.2 #gravitional force

def eq_of_motion(y,t,V0,gamma0):
    V1=y[0]
    gamma1=y[1]
    dV1dt=(2*g*sin(gamma0)/V0) * V1 - (g*cos(gamma0)) * gamma1
    dgamma1dt=(2*g*cos(gamma0)/(V0**2)) * V1 + (g*sin(gamma0)/V0) *gamma1
    
    return [dV1dt, dgamma1dt]#,dEdt,dhdt]


def yt(dydt,y0,dt):
    y=np.zeros(len(dydt))
    for i,dydti in enumerate(dydt):
        if i ==0:
            y[i]=y0
        else:
            y[i]=dydti*dt+y[i-1]
    return y


CL_LDmax=0.907
CD_LDmax=0.056
gamma0=np.arctan(-CD_LDmax/CL_LDmax)
V0=26
gammaInitial = [-0.1,0,0.1]
t=np.linspace(0,10,101)

fig,ax=plot.subplots(3,1,figsize=(10,20),gridspec_kw={'height_ratios': [3, 1,1]})


for gammai in gammaInitial:
    
    V0_gamma0=[0,gammai]

    y=odeint(eq_of_motion,V0_gamma0,t,args=(V0,gamma0,))
    
    V1=y[:,0]
    gamma1=y[:,1]
    V=V0+V1
    gamma=gamma0+gamma1
    
    dEdt=V*np.cos(gamma)
    dhdt1=V*np.sin(gamma)
    ht1=yt(dhdt1,300,t[1]-t[0])
    Et1=yt(dEdt,0,t[1]-t[0])
    ax[0].set_title('Intitial Flight Path condition')
    ax[0].plot(Et1,ht1)
    ax[0].set_ylabel('h (ft)')
    ax[0].set_xlabel(r'$\zeta$ (ft)')
    ax[1].plot(t,V)
    ax[1].set_ylabel('V (ft/s)')
    ax[1].set_xlabel(r'time (s)')
    ax[2].plot(t,gamma)
    ax[2].set_ylabel('gamma (deg)')
    ax[2].set_xlabel(r'time (s)')

ax[0].legend(['Downwards launch', 'Median launch', 'Upwards launch'])


fig,ax=plot.subplots(3,1,figsize=(10,20),gridspec_kw={'height_ratios': [3, 1,1]})
VInitial = [-1,0,1]


for Vi in VInitial:
    
    V0_gamma0=[Vi,gamma0]

    y=odeint(eq_of_motion,V0_gamma0,t,args=(V0,gamma0,))
    
    V1=y[:,0]
    gamma1=y[:,1]
    V=V0+V1
    gamma=gamma0+gamma1
    
    dEdt=V*np.cos(gamma)
    dhdt1=V*np.sin(gamma)
    ht1=yt(dhdt1,300,t[1]-t[0])
    Et1=yt(dEdt,0,t[1]-t[0])
    ax[0].set_title('Intitial Velocity condition')
    ax[0].plot(Et1,ht1)
    ax[0].set_ylabel('h (ft)')
    ax[0].set_xlabel(r'$\zeta$ (ft)')
    ax[1].plot(t,V)
    ax[1].set_ylabel('V (ft/s)')
    ax[1].set_xlabel(r'time (s)')
    ax[2].plot(t,gamma)
    ax[2].set_ylabel('gamma (deg)')
    ax[2].set_xlabel(r'time (s)')
    
ax[0].legend(['Low speed launch', 'Median speed launch', 'High speed launch'])