#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import tov_solver as tov
import matplotlib.pyplot as plt

def M_R_curve():
    
    """### PART A ###"""
    #### Scaling Factors ####
    length = 1477 #m
    time = 4.927 #s
    solar_mass = 1.989e30 #kg

    radius = []
    mass = []
    baryonic_mass = []

    rho_c = np.linspace(1e-1,1e-4,200) ## 200 evenly spaced rho_c values as inital guess around atomic density 1e-3

    ### PART A ###
    for rho in rho_c:
        sol = tov.TOV_solver(rho)
        radius.append(sol.t[-1])
        mass.append(sol.y[0,-1])
        baryonic_mass.append(sol.y[3,-1])

    mass = np.array(mass)
    radius = np.array(radius)
    baryonic_mass = np.array(baryonic_mass)
    
    radius = radius*(length/1000) #km
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(radius,mass,'g')
    plt.grid(1)
    plt.title("Mass-Radius of Neutron Stars")
    plt.xlabel("Radius(km)")
    plt.ylabel("Mass(Solar Mass)")
    rho_c = rho_c*(solar_mass/(length**3))
    
    """### PART B ###"""
    ## Fractional Binding energy calculation and plot ##
    frac_bind_energy = (baryonic_mass - mass)/mass
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(radius,frac_bind_energy,'g')
    plt.grid(1)
    plt.xlabel("Fractional Binding Energy")
    plt.ylabel("Radius(km)")
    
    """### PART C ###"""
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(rho_c,mass)
    plt.grid(1)
    plt.title("rho_c - mass")
    plt.xlabel("rho_c(kg/m3)")
    plt.ylabel("mass(solar_mass)")

    ###STABILTY###
    unstable_mass = []
    unstable_radius = []
    stable_mass = []
    stable_radius =  []
    tol = 1e-21
    for i in range(len(mass)):
        if (i == (len(mass)-1)):
            stable_mass.append(mass[i])
            stable_radius.append(radius[i])
        elif( ( (mass[i+1] - mass[i]) / (rho_c[i+1] - rho_c[i]) ) < 0):
            unstable_mass.append(mass[i])
            unstable_radius.append(radius[i])

        elif((mass[i+1] - mass[i]) / (rho_c[i+1] - rho_c[i])  <= tol):
            max_mass = mass[i]
            print("maximum mass is : ",max_mass)

        else:
            stable_mass.append(mass[i])
            stable_radius.append(radius[i])
                                     
    ### plot stability region ###
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(stable_radius,stable_mass,'--',unstable_radius,unstable_mass)
    plt.title("Stability of Neutron Stars")
    plt.xlabel("Radius(km)")
    plt.ylabel("Mass(Solar Mass)")
    plt.gca().legend(('stable NS',"unstable NS"), loc='upper right')
    plt.grid(1)

