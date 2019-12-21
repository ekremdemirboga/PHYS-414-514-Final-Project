#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import matplotlib.pyplot as plt

def chan_solve(D,rho_c):
    
    def f(x, y, z):
        return z

    def g(x, y, z, rho_c,D):
        y_c = np.sqrt((rho_c/D) + 1)
        return -(y**2 - 1/(y_c**2))**(3/2) -2/x*z

    def integrate(x_0, y_0, z_0, rho_c, D, t_step):
        '''Integrates one time step'''
        k_0 = t_step * f(x_0, y_0, z_0)
        l_0 = t_step * g(x_0, y_0, z_0, rho_c, D)
        k_1 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0)
        l_1 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0, rho_c, D)
        k_2 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1)
        l_2 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1, rho_c, D)
        k_3 = t_step * f(x_0+t_step, y_0+k_2, z_0+l_2)
        l_3 = t_step * g(x_0+t_step, y_0+k_2, z_0+l_2, rho_c, D)
        x_1 = x_0 + t_step
        y_1 = y_0 + 1/6 * (k_0+2*k_1+2*k_2+k_3)
        z_1 = z_0 + 1/6 * (l_0+2*l_1+2*l_2+l_3)
        return (x_1, y_1, z_1)

    # define initial conditions
    x_0 = 1e-50 # need some small value that is close to 0 but not exactly 0
    y_0 = 1
    z_0 = 0

    # do the integration and compile the results into two lists
    xs = [x_0]
    ys = [y_0]
    rho_c = 1e8
    D = 1e10
    t_step = 0.01
    y_c = np.sqrt((rho_c/D) + 1)
    while y_0>(1/((y_c)**2)):
        x_0, y_0, z_0 = integrate(x_0, y_0, z_0, rho_c, D, t_step)
        xs.append(x_0)
        ys.append(y_0)
    return xs,ys

