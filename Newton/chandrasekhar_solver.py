#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
from scipy.integrate import solve_ivp

def solve_chan(D,rho_c,limit=5,met = 'LSODA'):
    ##### takes the guess D and inital guess rho_c                                    ######
    ##### Then, solves the chandrasekhar white dwarf equation and returns it as solve.######
    ##### also returns the value which unction gives the surface for further use     ######
    def f(t, y):
        y_c = np.sqrt((rho_c/D) + 1)
        dydt = [y[1] , (-2/t)*y[1] - (y[0]**2 - (1/y_c**2))**(3/2)]
        return dydt
        
    y_c = np.sqrt((rho_c/D) + 1)
    dt = 0.01
    tspan = np.linspace(0+dt,limit,100)  
    yinit = [1,0]
    
    sol = solve_ivp(lambda t, y: f(t, y),[tspan[0], tspan[-1]], yinit,method=met,t_eval=tspan)
    surface = 1/(y_c)
    return sol,surface





