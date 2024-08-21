# -*- coding: utf-8 -*-
"""
Created on Sun May  5 13:05:14 2024

@author: tanne
"""
import numpy as np
import matplotlib.pyplot as plt

# Default Parameter values
r = 2.5
kappa = 68040
beta = 0.81
seed = 140000

# Base settings and labels
a_values = [2.0, np.e, 10.0]
a_values_str = ['2', 'e', '10']
colors = ['tab:blue', 'tab:orange', 'tab:green']

# Function to return the effective phi value based on base value (a)
def q_effective(a, phi):
    return np.log(phi + 1) / np.log(a)

# Define domain
phi = np.linspace(0.0, 0.9999, 50)

# Plot the functions
plt.figure(figsize=(4, 3))

# Loop over each value of a and plot the corresponding stable steady state value
for a, color, a_str in zip(a_values, colors, a_values_str):
    xi = q_effective(a, phi)
    y = [(beta * seed / kappa) * (1 - (1 / ((1 - xi_val) * r))) for xi_val in xi]
    
    # Replace any y values < 0 with 0 to only plot stable steady states
    y = [0 if val < 0 else val for val in y]
    
    plt.plot(phi, y, color=color, label=f'${a_str}$')

# Labeling plot and properties
plt.xlabel('$\phi$ - Fungicide Application Rate')
plt.ylabel('Stable Mold Steady State')
plt.legend(title='$a$', frameon=0, title_fontsize = '12')
plt.grid(False)

# Saving the plot
plt.savefig('bifurcation.png', dpi=600, bbox_inches="tight")

# Display the plot
plt.show()


