#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from math import pi, exp, sin
from scipy import fft


# In[2]:


def sinc(x):
    if x == 0:
        return 1
    else:
        return sin(x)/x


# In[3]:


def Source_Correctiion(original_time_array, original_dt, original_source_array, correction_time_array, correction_dt):
    # Forward Transform
    
    omega_range =  pi / original_dt
    omega_array = np.linspace(-omega_range, omega_range, original_time_array.size+1)
    domega = 2 * pi / (original_dt * original_time_array.size)
    
    aux_array = np.zeros(omega_array.size, dtype=complex)
    for index, item in enumerate(omega_array):
        tmp = np.exp(complex(0,1) * item * sinc(item * original_dt * 0.5) * original_time_array) * original_source_array * original_dt
        aux_array[index] = tmp.sum()
        
    # Inverse Transform
    correction_vibration_array = np.zeros(correction_time_array.size, dtype=complex)
    
    for index, item in enumerate(correction_time_array):
        tmp = np.exp(-1 * complex(0,1) * omega_array * item) * aux_array * domega
        correction_vibration_array[index] = tmp.sum()
        
    correction_vibration_array = correction_vibration_array.real * 0.5 / pi
    
    return correction_vibration_array


# In[4]:


def Wave_Correctiion(original_time_array, original_dt, original_source_array, correction_time_array, correction_dt):
    # Forward Transform
    
    omega_range =  pi / original_dt
    omega_array = np.linspace(-omega_range, omega_range, original_time_array.size+1)
    domega = 2 * pi / (original_dt * original_time_array.size)
    
    aux_array = np.zeros(omega_array.size, dtype=complex)
    for index, item in enumerate(omega_array):
        tmp = np.exp(complex(0,1) * item * sinc(item * original_dt * 0.5) * original_time_array) * original_source_array * original_dt
        aux_array[index] = tmp.sum()
        
    # Inverse Transform
    correction_vibration_array = np.zeros(correction_time_array.size, dtype=complex)
    
    for index, item in enumerate(correction_time_array):
        tmp = np.exp(-1 * complex(0,1) * omega_array * item) * aux_array * domega * np.cos(omega_array * original_dt * 0.5)
        correction_vibration_array[index] = tmp.sum()
        
    correction_vibration_array = correction_vibration_array.real * 0.5 / pi
    
    return correction_vibration_array


# In[ ]:




