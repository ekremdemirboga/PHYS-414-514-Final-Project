#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import tov_solver as tov
from scipy import interpolate


def wave_equation(rho_c):
    
    def M(x):
        ## interpolation of M(r) from TOV solution
        sol = tov.TOV_solver(rho_c)
        r = np.linspace(0,sol.t[-1],len(sol.y[0]))
        Mr = interpolate.interp1d(r, sol.y[0])
        return Mr(x)
    def v(x):
        ## interpolation of v(r) from TOV solution
        sol = tov.TOV_solver(rho_c)
        r = np.linspace(0,sol.t[-1],len(sol.y[0]))
        vr = interpolate.interp1d(r, sol.y[1])
        return vr(x)

    def func(r):
        ## function f(r) in the wave equation
        sol = tov.TOV_solver(rho_c)
        return np.exp(M(r)/2)*(1 - 2*v(r)/r)**(1/2)

    def gfunc(r):
        ## defined function g for wave equation
        return 1/r**2

    def init(r):
        ## inital condition for the wave
        return -2*1e-3*np.exp(-r**2/2)*r

    sol = tov.TOV_solver(rho_c)
    c = 1 #speed of light
    N = 300 # mesh element number
    r = np.linspace(0,sol.t[-1],N) #space
    t = np.linspace(0,5,N) #time
    dr = r[1] - r[0] #differential space element
    dt = t[1] - t[0] # differential time element
    lam = c*dt/dr # lambda

    #discretize functions
    g = np.zeros(N)
    g = gfunc(r)
    f = np.zeros(N)
    f = func(r) 
    f[0] = 1
    g[0] = 1e5
    Psi = np.zeros((N,N))
    Pi = np.zeros((N,N))

    #initial condition
    for i in range(N):
        Psi[0,i] = init(r[i])
        Pi[0,i] = init(r[i])/(1/f[i])

    #BCs
    Psi[:,0] = 0
    Pi [:,0] = 0

    Psi[:,-1] = 0
    Pi [:,-1] = 0

    #solving PDE
    T = 200
    for n in range(1,T): ##time
        for j in range(1,N-1): ##space
            if (n == 0):
                Psi[1,j] = 0.75*Psi[0,j] + 0.25*Psi[2,j] + lam*0.5*((f*Pi)[0,j+1] - (f*Pi)[0,j-1])
                Pi[1,j] = 0.75*Pi[0,j] + 0.25*Pi[2,j] + lam*0.5*g[j]*(((1/g)*f*Psi)[0,j+1] - ((1/g)*f*Psi)[0,j-1])

            Psi[n+1,j] = Psi[n-1,j] + lam*((f[j+1]*Pi[n,j+1]) - (f[j+1]*Pi[n,j-1]))
            Pi[n+1,j] = Pi[n-1,j] + lam*g[j]*(((1/g)*f*Psi)[n,j+1] - ((1/g)*f*Psi)[n,j-1])         
    #plotting space vs time

    
    psi = np.zeros((N,N))
    psi[:,0] = 1e-3
    for l in range(N-1):
        for n in range(T):
            psi[n,l+1] = psi[n,l] + dr*Psi[n,l+1]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(r[2:], t[0:T:20])
    Z = psi[0:T:20,2:]
    ax.set_xlabel('r')
    ax.set_ylabel('time')
    ax.set_zlabel('Psi')
    surf = ax.scatter3D(X, Y, Z, cmap='viridis')
    ax.set_title('wave eq');       

