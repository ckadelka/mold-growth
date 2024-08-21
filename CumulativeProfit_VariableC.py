# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:26:46 2024

@author: tanne
"""
import numpy as np
import matplotlib.pyplot as plt

# Defualt parameters
x0 = 0.1
r = 2.5
kappa = 68040
beta = 0.81
T = 10
d = 1 / 1744
b = 12
theta = 0.7
seed = 140000

# Define different c values for each a
c_values_list = [
    [0.00547, 0.00647, 0.00747],  # c values for a = 2
    [0.00545, 0.00645, 0.00745],  # c values for a = e
    [0.0047, 0.0051, 0.0055]   # c values for a = 10
]

phi_values = np.linspace(0.0, 1, 1000)  # Avoid zero phi

# Define colors for different a values and line styles for c values
colors = ['tab:blue', 'tab:orange', 'tab:green']
line_styles = ['-', '--', ':']

def q_effective(a, phi):
    return np.log(phi + 1) / np.log(a)

def get_mstar(beta, s, kappa, phi, r):
    return (beta * s / kappa) * (1 - (1 / ((1 - phi) * r)))

def get_profit(phi, r, kappa, beta, seed, c, theta, m, b):
    m_0 = get_mstar(beta, seed, kappa, 0, r)
    profit = b * (d * beta * seed) * (1 - theta * (m / m_0)) - c * phi * beta * seed
    return profit

def get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, c):
    x = [x0]
    x_old = x0
    profit = 0

    for i in range(T):
        x_new = (1 - phi) * r * x_old * (1 - kappa * (x_old / (beta * seed)))
        x.append(x_new)
        profit += get_profit(phi, r, kappa, beta, seed, c, theta, x_new, b)
        x_old = x_new
    return profit

# Subplots for different a values
f, ax = plt.subplots(1, 3, figsize=(11, 4.5))
a_values_str = ['2', 'e', '10']

#plot all lines for requisite conditions
for ii, (a, c_values) in enumerate(zip([2, np.e, 10], c_values_list)):
    for c_idx, c in enumerate(c_values):
        profits = []
        for phi_1 in phi_values:
            phi = q_effective(a, phi_1)
            profit = get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, c)
            profits.append(profit)
        ax[ii].plot(phi_values, profits, label=f'c = {c:.5f}', linestyle=line_styles[c_idx], color=colors[ii])

    #plot line at profit when applying no fungicide S
    profit_at_phi_0 = get_timeseries_profit(x0, 0, r, kappa, beta, seed, T, c_values[0])
    ax[ii].axhline(y=profit_at_phi_0, linestyle='--', color='black', linewidth=1)
    ax[ii].legend(title='Fungicide Cost ($/plant)', frameon=0, loc='best')
    if ii == 0:
        ax[ii].set_ylabel('Cumulative Profit per Acre ($)')
    ax[ii].set_xlabel('$\phi$ - Fungicide Application Rate')

    # Make all spines (box) visible
    ax[ii].spines['right'].set_visible(True)
    ax[ii].spines['top'].set_visible(True)
    ax[ii].spines['left'].set_visible(True)
    ax[ii].spines['bottom'].set_visible(True)

plt.subplots_adjust(wspace=0.3)
plt.savefig('explanatory_plot_updated.pdf', bbox_inches='tight')
plt.show()



