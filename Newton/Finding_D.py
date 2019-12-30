#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import chandrasekhar_solver as ch
import Functions as F
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def Find_D():
    ### DATA ### 
    M,logg = F.data_reader()
    M,R = F.scaler(M,logg)

    ### inital guess ###
    max_iter = 200
    D = 1700000000.0
    dD = 1e7
    rho_c = np.linspace(10**10,10**7,20)
    tol = 7e-9

    ### CONSTANTS(SI) ###
    K = 3144530.473379261 #found before
    G = 6.67408e-11
    solar_mass = 1.988e30 #kg
    earth_radius = 6.371e6 #m

    zeta_n = 0 ## when the solution reaches surface value.
    zeta_prime_n = 0 ## when the derivative reaches the surface
    Radius = 0
    Mass = 0
    R_new = [] ## for storing new radius and mass
    M_new = []


    for i in range(max_iter):

        #D or rho_c dependent variables
        C = (5/8)*K*(D**(5/3))
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
                


            #############################
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
        R_new = np.array(R_new)
        M_new = np.array(M_new)
        spl = UnivariateSpline(R_new, M_new)
        R_new = []
        M_new = []
    
        ### calculate error
        error=0
        for l in range(len(R)):
            error = error + (M[l]-spl(R[l]))**2
            error = error/(len(R))
    
        if (error < tol):
            print("D is :",D,"with error: ", error,"corresponding C is: ", C)
            x = np.linspace(0,2.8,500)
            plt.figure()
            fig, axis = plt.subplots(figsize = (9,5))
            plt.plot(R,M,'ro',x,spl(x))
            plt.title("Interpolation of mass-radius")
            plt.xlabel("Radius(earth radius)")
            plt.ylabel("Mass(solar mass)")
            plt.gca().legend(('Data','interpolation for given D'), loc='upper right')
            break;
      
        """CHECK FIT look if D is okay, if not update it """
        if (sum(M)-sum(spl(R))>0):
            D = D - dD
        else:
            D = D + dD      
        print(i,"error:",error," and corresponding D: ",D)
        error_ref = error
        

