import csv
import numpy as np
import matplotlib.pyplot as plt

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
                
    star_masses = np.array(star_masses)
    star_surface_g = np.array(star_surface_g)
    

    return star_masses,star_surface_g

def scaler(M,logg):
    ####takes mass in solar mass and log of g, return scaled mass and radius ####
    ####in solar mass and earth radii####
    #constants#
    solar_mass = 1.988e30 #kg
    G = 6.67408e-11
    earth_radius = 6.371e6 #m

    g = (10**logg)*(10**(-2)) # in SI units
    mass = M * solar_mass # mass in kg
    R = np.sqrt(G*mass/g) #m
    R = R/earth_radius
    return M,R

def low_mass_radius(M,R,max_mass):
    #### takes arrays M and R returns the mass array M_new with smaller than max_mass ####
    #### and corresponding R_new arrays ####
    M_new = []
    R_new = []
    for i in range(len(M)):
        if (M[i] <= max_mass):
            M_new.append(M[i])
            R_new.append(R[i])
    R_new = np.array(R_new)
    M_new = np.array(M_new)
    return M_new,R_new
        
