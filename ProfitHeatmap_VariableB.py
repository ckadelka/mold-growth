# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 12:48:25 2024

@author: tanne
"""

import numpy as np
import matplotlib.pyplot as plt


#  Default arameters
x0 = 0.1
r = 2.5
kappa = 68040
beta = 0.81
T = 10
#plant to bushel conversion
c = 1 / 1744
theta = 0.7
d = 0.002
seed = 140000

phi_values = np.linspace(0, 1, 100)
b_values = np.linspace(0, 20, 31)
a_values = [2, np.e, 10]
a_str = ['2', 'e', '10']


# Custom color maps and marker colors
colors = {
    2: 'Blues_r',
    np.e: 'Oranges_r',
    10: 'Greens_r'
}
marker_colors = {
    2: 'black',
    np.e: 'black',
    10: 'black'
}

def q_effective(a, phi):
    return np.log(phi + 1) / np.log(a)

#calculate steady state of mold 
def get_mstar(beta, s, kappa, phi, r):
    return (beta * s / kappa) * (1 - (1 / ((1 - phi) * r)))

#calculate profit for a single year
def get_profit(phi, r, kappa, beta, seed, c, theta, m, b):
    #calculate how bad the mold could get
    m_0 = get_mstar(beta, seed, kappa, 0, r)
    profit = b * (c * beta * seed) * (1 - theta * (m / m_0)) - d * phi * beta * seed
    return profit

def get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, b):
    x = [x0]
    x_old = x0
    profit = 0

    for i in range(T):
        x_new = (1 - phi) * r * x_old * (1 - kappa * (x_old / (beta * seed)))
        x.append(x_new)
        profit += get_profit(phi, r, kappa, beta, seed, c, theta, x_new, b)
        x_old = x_new

    return profit

fig, axs = plt.subplots(3, 1, figsize=(5, 10))

for idx, a in enumerate(a_values):
    # Initialize an empty matrix to store profit values
    p_matrix = np.zeros((len(phi_values), len(b_values)))

    # Calculate profit for each combination of phi and b
    for i, phi_1 in enumerate(phi_values):
        phi = q_effective(a, phi_1)
        for j, b in enumerate(b_values):
            p = get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, b)
            p_matrix[i, j] = p

    ax = axs[idx]

    # Create the heatmap with custom color maps
    im = ax.imshow(p_matrix, cmap=colors[a], interpolation='nearest', origin='lower', aspect='auto')
    fig.colorbar(im, ax=ax, label='Cumulative Profit per Acre ($)')

    # Set the tick positions and labels for the x-axis (b values)
    xtick_positions = np.linspace(0, len(b_values) - 1, num=5, dtype=int)
    xtick_labels = np.around(b_values[xtick_positions], 0)
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels)

    # Set the tick positions and labels for the y-axis (phi values)
    ytick_positions = np.linspace(0, len(phi_values) - 1, num=6, dtype=int)
    ytick_labels = np.around(phi_values[ytick_positions], 1)
    ax.set_yticks(ytick_positions)
    ax.set_yticklabels(ytick_labels)

    # Hide x-axis labels for the first two subplots
    if idx < len(a_values) - 1:
        ax.set_xticklabels([])

    # Find the optimal phi for each b value and plot the markers
    optimal_phi_indices = np.argmax(p_matrix, axis=0)
    for j in range(len(b_values)):
        ax.scatter(j, optimal_phi_indices[j], color=marker_colors[a], marker='.', s=60)

    # Add labels and title
    ax.set_ylabel('$\phi$ - Fungicide Application Rate')
    legend_label = 'Optimal $\phi$'
    ax.legend([plt.Line2D([0], [0], color=marker_colors[a], marker='.', linestyle='', markersize=10)], [legend_label], loc='lower left', bbox_to_anchor=(0.55, 0.0))

# Only show x-axis label for the bottom subplot
axs[-1].set_xlabel('b - Bushel Price ($)')


# Create an additional axis for the label
fig.subplots_adjust(left=0.1, right=0.85)
ax_left = fig.add_axes([-0.1, 0.12, 0.02, 0.76])
ax_left.set_xticks([])
ax_left.set_yticks([])
for spine in ax_left.spines.values():
    spine.set_visible(False)
ax_left.plot([0, 0], [0, 1], 'k-', lw=0.5)
ax_left.set_ylim([0, 1])
ax_left.text(-0.15,0.5,'Maximal Fungicide Efficiency',ha='center',va='center',rotation=90)  
ax_left.text(0.15, 0.166, 'Low ($a=10$)', ha='center', va='center', rotation=90)
ax_left.text(0.15, 0.5, 'Medium ($a=e$)', ha='center', va='center', rotation=90)
ax_left.text(0.15, 0.833, 'High ($a=2$)', ha='center', va='center', rotation=90)

fig.subplots_adjust(wspace=0.05, hspace=0.05)
fig.savefig('fig2_profitheatmap_bushel.pdf',bbox_inches = "tight") 
plt.show()
