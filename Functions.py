#!/usr/bin/env python
# coding: utf-8

# In[2]:


import csv
import numpy as np

def data_reader(file_name='white_dwarf_data.csv'):
    ### Reads the file and returns columns as sorted arrays ###
    star_masses = []
    star_surface_g = []
    line_count = 0
    with open(file_name) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            if line_count == 0:
                line_count +=1
            else:
                star_surface_g.append(float(row[1]))
                star_masses.append(float(row[2]))
                line_count += 1
                
    star_masses = np.sort(np.array(star_masses))
    star_surface_g = np.sort(np.array(star_surface_g))
    

    return star_id,star_masses,star_surface_g
            
def plot_m_vs_R(M,R):
    #constants
    solar_mass = 1.988e30 #kg
    G = 6.67408e-11
    earth_radius = 6.371e6 #m
    
    star,mass,logg = data_reader()
    g = (10**logg)*(10**(-2)) # in SI units
    mass = mass * solar_mass # kg
    R = np.sqrt(G*mass/g) #m
    
    plt.xlabel("Mass of the stars(solar_mass)")
    plt.ylabel("Radius of the stars(earth radii)")
    plt.title("Mass vs Radius for WDs")
    
    figure = plt.plot(mass/solar_mass,R/earth_radius,'o')
    return figure


# In[ ]:




