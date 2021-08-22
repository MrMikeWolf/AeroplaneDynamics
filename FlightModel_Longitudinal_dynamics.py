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


"""
BOEING 747
"""
q0bar=177
S=5500
m=19790
V0=774
c=27.31
CLalpha=4.9
CLalphadot=-5.978
CL0=0.65
CLq=5.965
CLd=0.37
g=32.2# ft/s
Cmadot=-6.4
Iyy=3.331e7
Cma=-1.03
Cmq=-24
Cmd=-1.44
CDalpha=0.43
CD= lambda alpha : alpha*0.0000001 # aerodynamic coefficient of drag : f(alpha)
rho=1 # i think density of the air

#UNKNOWN
Zc=1
Mc=1

Za_dot=-q0bar*S*c*CLalphadot/(2*m*V0**2)

Za=-q0bar*S*CLalpha/(m*V0)
Zq=-q0bar*S*c*CLq/(2*m*V0**2)
Zd=-q0bar*S*CLd/(m*V0)

Ma_dot=q0bar*S*c**2*Cmadot/(2*Iyy*V0)

Ma=q0bar*S*c*Cma/Iyy
Mq=q0bar*S*c**2*Cmq/(2*Iyy*V0)
Md=q0bar*S*c*Cmd/Iyy

Za_apost=Za/(1-Za_dot)
Zq_apost=(Zq+Za_dot)/(1-Za_dot)
Zd_apost=(Zd)/(1-Za_dot)

Ma_apost=Ma+Ma_dot*Za/(1-Za_dot)
Mq_apost=Mq+Ma_dot*(1+Zq)/(1-Za_dot)
Md_apost=Md+Ma_dot*Zd/(1-Za_dot)

alpha0=(-1/(Za*Md-Ma*Zd))*((g/V0+Zc)*Md-Mc*Zd)
de0=(-1/(Za*Md-Ma*Zd))*(Za*Mc-Ma*(g/V0+Zc))

#Overwrite values as Zc and Mc are not given
alpha0=4.6
de0=(m*g/(q0bar*S)-CL0-CLalpha*alpha0)

#initial
T0=q0bar*S*CD(alpha0)/np.cos(alpha0)
gamma0=-0.1

def eq_of_motion(y,t,V0,gamma0,alpha0,de0):
    V1=y[0]
    gamma1=y[1]
    alpha1=y[2]
    q1=y[3]
    
    dV1dt=(-rho*V0*S*CD(alpha0))/m * V1 - g * gamma1
    dgamma1dt=(2*g/(V0**2)-2*T0*sin(alpha0)/(m*V0**2)) * V1 + 0 * gamma1 # REMOVE THE 15 (user adjusted)
    dalpha1dt=Za_apost*alpha1+ (1+Zq_apost)*q1 +Zd_apost*de0
    dq1dt=Ma_apost*alpha1+Mq_apost*q1+Md_apost*de0
    
    return [dV1dt, dgamma1dt, dalpha1dt, dq1dt]

def yt(dydt,y0,dt):
    y=np.zeros(len(dydt))
    for i,dydti in enumerate(dydt):
        if i ==0:
            y[i]=y0
        else:
            y[i]=dydti*dt+y[i-1]
    return y

# # INITIAL CONDITIONS
t=np.linspace(0,200,101)

fig,ax=plot.subplots(7,1,figsize=(10,10))#,gridspec_kw={'height_ratios': [3, 1,1]})

V0_gamma0_alpha0_q0=[0,gamma0,alpha0,q0bar]

y=odeint(eq_of_motion,V0_gamma0_alpha0_q0,t,args=(V0,gamma0,alpha0,de0))

V1=y[:,0]
gamma1=y[:,1]
V=V0+V1
gamma=gamma0+gamma1
dEdt=V*np.cos(gamma)
dhdt1=V*np.sin(gamma)
ht1=yt(dhdt1,100,t[1]-t[0])
Et1=yt(dEdt,0,t[1]-t[0])
theta=y[:,2]+y[:,1]

ax[0].set_title('Intitial Velocity condition')
ax[0].set_ylabel('V (ft/s)')
ax[0].plot(t,y[:,0])
ax[1].set_ylabel(r'$\gamma$ (deg)')
ax[1].plot(t,y[:,1])
ax[2].set_ylabel(r'$\alpha$ (deg)')
ax[2].plot(t,y[:,2])
ax[3].set_ylabel('q (deg/s)')
ax[3].plot(t,y[:,3])
ax[4].set_ylabel(r'$\theta$ (deg)')
ax[4].plot(t,theta)
ax[5].set_ylabel(r'h (ft)')
ax[5].plot(t,ht1)
ax[6].set_ylabel(r'$\zeta$ (ft)')
ax[6].plot(t,Et1)
fig.tight_layout()
plot.show()
