#!/usr/bin/env python
# coding: utf-8

# In[23]:


import numpy as np
### THEORATİCAL VALUES OF C AND D###

hbar = 1.0545718e-34
c = 299792458
me = 9.10938356e-31
atomic_mass_unit = 1.66054e-27
G = 6.674081e-11
solar_mass = 1.988e30

C = ((me**4)*(c**5))/(24*(np.pi**2)*(hbar**3))
D = ((3*(np.pi**2)*(hbar**3))/(atomic_mass_unit*(me**3)*(c**3)*2))**(-1)


# In[24]:


print("Theoratical value of C is :",C)


# In[25]:


print("Theoratical value of D is :",D)


# In[26]:


K = (8*C)/(5*(D**(5/3)))


# In[27]:


print("Theoratical value of K is :",K)


# In[28]:


### Theorica value of Chandrasekhar
xi_3 = 6.926992292292292
theta_3_prime=-0.04210902023243969
K1 = (2*C/D**(4/3))
Mch = 4*np.pi*(K1/(G*np.pi))**(3/2)*xi_3**2*(-theta_3_prime)


# In[29]:


print("Theoratical value of Chandrasekhar mass is :",Mch/solar_mass)

