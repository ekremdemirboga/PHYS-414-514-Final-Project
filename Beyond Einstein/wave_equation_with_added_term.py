#!/usr/bin/env python
# coding: utf-8

# In[24]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import tov_solver as tov
from scipy import interpolate

def modified_wave_equation(rho_c = 1e-3):
    
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
    def p(x):
        ## interpolation of p(r) from TOV solution
        sol = tov.TOV_solver(rho_c)
        r = np.linspace(0,sol.t[-1],len(sol.y[0]))
        pr = interpolate.interp1d(r, sol.y[2])
        return pr(x)

    def rho(x):
        ## interpolation of rho(r) from TOV solution
        return (p(x)/50)**(1/2)

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
    def init2(r):
        ## inital condition for the wave
        return 1e-3*np.exp(-r**2/2)

    def hfunc(r):
        return rho(r)-3*p(r)

    sol = tov.TOV_solver(rho_c)
    c = 1 #speed of light
    N = 300 # mesh element number
    r = np.linspace(0,7,N) #space
    t = np.linspace(0,5,N) #time
    dr = r[1] - r[0] #differential space element
    dt = t[1] - t[0] # differential time element
    lam = c*dt/dr # lambda

    #discretize functions
    g = np.zeros(N)
    g = gfunc(r)
    f = np.zeros(N)
    f = func(r)
    h = np.zeros(N)
    h = hfunc(r)
    f[0] = 1
    g[0] = 1e4
    Psi = np.zeros((N,N))
    Pi = np.zeros((N,N))
    psi = np.zeros((N,N))
    #initial condition
    for i in range(N):
        Psi[0,i] = init(r[i])
        psi[0,i] = init2(r[i])

    #BCs
    Psi[:,0] = 0
    Pi [:,0] = 0
    psi[:,0] = 1e-3

    Psi[:,-1] = 0
    Pi [:,-1] = 0

    #solving PDE
    T = 250
    for n in range(1,T): ##time
        for j in range(1,N-1): ##space
            if (n == 0):
                Psi[1,j] = 0.75*Psi[0,j] + 0.25*Psi[2,j] + lam*0.5*((f[j+1]*Pi[0,j+1]) - (f[j-1]*Pi[0,j-1]))
                Pi[1,j] = 0.75*Pi[0,j] + 0.25*Pi[2,j] + lam*0.5*g[j]*(((1/g[j+1])*f[j+1]*Psi[0,j+1]) - ((1/g[j-1])*f[j-1]*Psi[0,j-1]))
                +48*np.pi*np.exp(-12*psi[0,j]**2)*psi[0,j]*(h[j])*dt
                psi[n,j+1] = psi[n,j] + dr*Psi[n,j+1]

            Psi[n+1,j] = Psi[n-1,j] + lam*((f[j+1]*Pi[n,j+1]) - (f[j+1]*Pi[n,j-1]))
            Pi[n+1,j] = Pi[n-1,j] + lam*g[j]*(((1/g[j+1])*f[j+1]*Psi[n,j+1]) - ((1/g[j-1])*f[j-1]*Psi[n,j-1])) 
            +48*np.pi*np.exp(-12*psi[n,j]**2)*psi[n,j]*(h[j])*dt
            psi[n,j+1] = psi[n,j] + dr*Psi[n,j+1]
    #plotting space vs time
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(r[2:], t[0:T:20])
    Z = psi[0:T:20,2:]
    ax.set_xlabel('r')
    ax.set_ylabel('time')
    ax.set_zlabel('psi')
    surf = ax.scatter3D(X, Y, Z, cmap='viridis')
    ax.set_title('Solution for psi');
    ax.set_xlim3d(0, sol.t[-1])
    #ax.set_ylim3d(0,0.2)
    ax.set_zlim3d(-0.001,0.000)

