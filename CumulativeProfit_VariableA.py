# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:21:53 2024

@author: tanne
"""
import numpy as np
import matplotlib.pyplot as plt

# Default Parameters
x0 = 0.1
r = 2.5
kappa = 68040
beta = 0.81
T = 10
#conversion rate(plant to bundle) x $12 bundle price
c = 12 / 1744
theta = 0.7
d = 0
seed = 140000

# Plot settings
a_values = [2.0, np.e, 10.0]
a_values_str = ['$2$', '$e$', '$10$']
colors = ['tab:blue','tab:orange','tab:green']
phi_values = np.linspace(0, 1, 100)

# Function Definitions
def q_effective(a, phi):
    return np.log(phi + 1) / np.log(a)

def get_mstar(beta, s, kappa, phi, r):
    return (beta * s / kappa) * (1 - (1 / ((1 - phi) * r)))

def get_profit(phi, r, kappa, beta, seed, c, theta, m, d):
    m_0 = get_mstar(beta, seed, kappa, 0, r)
    profit = (c * beta * seed) * (1 - theta * (m / m_0)) - d * phi * beta * seed
    return profit

def get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, d):
    x = [x0]
    x_old = x0
    #total_profit = get_profit(phi, r, kappa, beta, seed, c, theta, x0, d)
    total_profit = 0

    for i in range(T):
        x_new = (1 - phi) * r * x_old * (1 - kappa * (x_old / (beta * seed)))
        x.append(x_new)
        total_profit += get_profit(phi, r, kappa, beta, seed, c, theta, x_new, d)
        x_old = x_new

    return total_profit

# Initialize the plot
fig, ax = plt.subplots(figsize=(4, 3))

# Calculate and plot profits for each value of a
for a, label, color in zip(a_values, a_values_str, colors):
    profits = []
    for phi_1 in phi_values:
        phi = q_effective(a, phi_1)
        profit = get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, d)
        profits.append(profit)
    plt.plot(phi_values, profits, label=label, color=color)

# Set plot labels and title
plt.xlabel('$\phi$ - Fungicide Application Rate')
plt.ylabel('Cumulative Revenue per Acre ($)')
#plt.title('Cumulative Profit for Different Values of $a$')
plt.legend(title='$a$', title_fontsize = '12', frameon=0)
ax.grid(False)
#ax.set_ylim(2500, 8000)

plt.savefig('fig1c.png', dpi=200, bbox_inches="tight")

# Show the plot
plt.show()
