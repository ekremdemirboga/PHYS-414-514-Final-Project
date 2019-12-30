#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
import chandrasekhar_solver as ch
import Functions as F
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def chandrasekhar_mass():
    
    ### DATA ### 
    M,logg = F.data_reader()
    M,R = F.scaler(M,logg)

    D = 1830000000 #found
    
    rho_c = np.linspace(1e20,1e22)
    ### CONSTANTS(SI) ###
    K = 3097671.1345466143 #found before
    G = 6.67408e-11
    solar_mass = 1.988e30 #kg
    earth_radius = 6.371e6 #m

    zeta_n = 0 ## when the solution reaches surface value.
    zeta_prime_n = 0 ## when the derivative reaches the surface
    Radius = 0
    Mass = 0
    R_new = [] ## for storing new radius and mass
    M_new = []
    convergence = []
    C = 5.165215436316856e+21 #found
    counter = 0
    for rho in rho_c:
        ### solve the chandrasekhar for given rho and D ###
        sol,surface = ch.solve_chan(D,rho,limit=8.2,met='RK45')
        ##### find the surface #####
        for j in range(len(sol.t)):
            if (sol.y[0,j] <= surface):
                zeta_n = sol.t[j]
                zeta_prime_n = sol.y[1,j]
            elif(sol.y[0,-1]> surface):
                zeta_n = sol.t[-1]
                zeta_prime_n = sol.y[1,-1]
        ### For Radius ###
        y_c = np.sqrt(rho/D + 1)
        Beta = np.sqrt((2*C)/(np.pi*G))/(D*y_c) ## scale factor of radius
        Radius = (Beta*zeta_n)
        ### For Mass ###
        Mass = 4*np.pi*(Radius**3)*D*(y_c**3)*(-zeta_prime_n/zeta_n)
        Mass = Mass/solar_mass
        Radius = Radius/earth_radius
        ## add them to array
        R_new.append(Radius)
        M_new.append(Mass)
        if (counter >0):
            convergence.append(M_new[-1]-M_new[-2])
        counter =+ 1
    R_new = np.array(R_new)
    M_new = np.array(M_new)
    convergence = np.array(convergence)
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(R_new,M_new,'r--',R_new,M_new,'o')
    plt.title("mass-radius")
    plt.xlabel("Radius(earth radius)")
    plt.ylabel("Mass(solar mass)")
    print("given rho_c's : ", rho_c)
    print("and corresponding masses: ",M_new)
    print("difference between two last mass values: ", M_new[-1]-M_new[-2])
    print("Chandrasekhar Mass is : ",M_new[-1],"(in solar mass)")
    plt.figure()
    fig, axis = plt.subplots(figsize = (9,5))
    plt.title("Convergence")
    plt.plot(convergence,'o')

