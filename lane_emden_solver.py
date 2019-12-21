#!/usr/bin/env python
# coding: utf-8

# In[100]:


import numpy as np
from scipy.integrate import solve_ivp

def solve_le(n,limit):
    ### Solves the Lane-Emden Equation for given n starting from 0 to limit. Uses 'RK45'.###
    ### Returns the result as scipy.integrate._ivp.ivp.OdeResult ###
    def f(t, y):
        a= n
        dydt = [y[1] , (-2/t)*y[1]-y[0]**(a)]  
        return dydt
    
    dt = 0.0001
    tspan = np.linspace(0+dt,limit,1000)  
    yinit = [1,0]
    
    sol = solve_ivp(lambda t, y: f(t, y),[tspan[0], tspan[-1]], yinit, t_eval=tspan)

    return sol

