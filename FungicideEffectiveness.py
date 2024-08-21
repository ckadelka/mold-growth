# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:29:31 2024

@author: tanne
"""
import numpy as np
import matplotlib.pyplot as plt

# Define the domain
phi = np.linspace(0, 1, 500)

# Define the three bases with colors and labels
a_values = [2, np.e, 10]
labels = [r'$a = 2$', r'$a = e$', r'$a = 10$']
colors = ['tab:blue','tab:orange','tab:green']

# Plot the functions
plt.figure(figsize=(4, 3))

for a, label, c in zip(a_values, labels, colors):
    
    #multiple by 100 to get percent mold reduction
    y = np.log(phi + 1) / np.log(a) * 100
    plt.plot(phi, y, label=label, color=c)

# Add plot properties
plt.xlabel(r'$\phi$ - Fungicide Application Rate')
plt.ylabel('$q(\phi)$ - Mold Reduction (%)') 
plt.legend(fontsize = '10',title=r'$q(\phi) = \log_a(\phi + 1)$', title_fontsize='12',frameon=0)
plt.grid(False)

# Save Plot
plt.savefig('fig1a.pdf',bbox_inches="tight")

# Show the plot
plt.show()

