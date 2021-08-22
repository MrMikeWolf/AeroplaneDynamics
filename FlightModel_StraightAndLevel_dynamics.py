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

def eq_of_motion(y,t,V0,gamma0):
    V1=y[0]
    gamma1=y[1]
    dV1dt=(-rho*V0*S*CD(alpha0))/m * V1 - g * gamma1
    dgamma1dt=(2*g/(V0**2)-2*T0*sin(alpha0)/(m*V0**2)) * V1 + 0 * gamma1
    
    return [dV1dt, dgamma1dt]#,dEdt,dhdt]

def yt(dydt,y0,dt):
    y=np.zeros(len(dydt))
    for i,dydti in enumerate(dydt):
        if i ==0:
            y[i]=y0
        else:
            y[i]=dydti*dt+y[i-1]
    return y

# INITIAL CONDITIONS
S = 1 # wing reference area
CL= lambda alpha : alpha-1  # aerodynamic coefficient of lift : f(alpha)
CD= lambda alpha : alpha*0.05 # aerodynamic coefficient of drag : f(alpha)
alpha0=0.05
m = 1 # mass of the aircraft (slugs)
g = 32.2# gravitional force
rho=1 # density of the air
# INITIAL CONDITIONS
V0=20
gamma0=0 #flight path angle is zero for straight and level Equilibirum

#UNKNOWNS: alpha and thrust
q0bar=0.5*rho*V0**2
zero=alpha0-np.arctan((m*g-q0bar*S*CL(alpha0))/q0bar*S*CD(alpha0))
T0=q0bar*S*CD(alpha0)/np.cos(alpha0)
xi=rho*V0**2*S*CD(alpha0)/(2*np.sqrt(2)*m*g)
t=np.linspace(0,20,101)
gammaInitial = [-0.1,0,0.1]

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
    ht1=yt(dhdt1,100,t[1]-t[0])
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