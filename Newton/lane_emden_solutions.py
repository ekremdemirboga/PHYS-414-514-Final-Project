#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from lane_emden_solver import solve_le


# In[3]:


def solutions15():    
    ### For n=1.5 ###
    fig, axis = plt.subplots(figsize = (9,5))
    sol = solve_le(3/2,5)
    plt.plot(sol.t,sol.y[0])
    plt.title("solution of Lane-Emden Equation for dn=1.5")
    plt.xlabel("Scaled Radius")
    plt.ylabel("Dimensionless Density")
    xi_n = sol.t[-1]
    print("xi_n for n=1.5",xi_n)
    plt.figure()
    fig, axis = plt.subplots(figsize = (9,5))
    plt.plot(sol.t,sol.y[1],'r')
    plt.title("derivative of solution of Lane-Emden Equation for dn=1.5")
    plt.xlabel("Scaled Radius")
    plt.ylabel("Theta'")
    theta_prime_n = sol.y[1][-1]
    print("theta'(xi_n) is ",theta_prime_n)


# In[2]:


def solutions():
    #for different n values
    n= [1,2,3,4]
    fig, axis = plt.subplots(figsize = (9,5))
    for i in n:
        sol = solve_le(i,20)
        plt.plot(sol.t,sol.y[0])
        print(i)
        for j in range(len(sol.t)):
            if (sol.y[0,j]<= 0):
                print("xi_n for n={}:".format(i), sol.t[j])
                print("xi_n prime for n={}:".format(i),sol.y[1,j])
                break
    plt.ylim((0,1))
    plt.xlim((0,20))
    plt.gca().legend(('n=1','n=2','n=3','n=4','n=5'), loc='upper right')
    plt.title("solutions of Lane-Emden Equation for different n values")
    plt.xlabel("Scaled Radius")
    plt.ylabel("Dimensionless Density")

