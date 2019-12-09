#!/usr/bin/env python
# coding: utf-8

# In[21]:


import csv
import numpy as np

def data_reader(file_name='white_dwarf_data.csv'):
    ### Reads the file and returns columns as arrays ###
    star_id = []
    star_masses = []
    star_surface_g = []
    line_count = 0
    with open(file_name) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            if line_count == 0:
                line_count +=1
            else:
                star_id.append(row[0])
                star_surface_g.append(float(row[1]))
                star_masses.append(float(row[2]))
                line_count += 1
                
    star_id = np.array(star_id)
    star_masses = np.array(star_masses)
    star_surface_g = np.array(star_surface_g)

    return star_id,star_masses,star_surface_g
            
            


# In[ ]:




