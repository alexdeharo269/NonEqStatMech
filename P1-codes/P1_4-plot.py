# Represents 2 plots with first and last time stamps of the brownian particles
import seaborn as sns
import numpy as np  
import matplotlib.pyplot as plt

# We read from a file which has the x and y positions of all particles at each time step
data = np.loadtxt('P1_4-data.txt')
# Set seaborn style for professional appearance
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)
sns.set_palette("husl")
#We read all rows since there are only 2 time steps, and all columns since we have all particles
x_start = data[0, ::2]  # x positions at first time step (even indices)
y_start = data[0, 1::2] # y positions at first time step (odd indices)
x_end = data[-1, ::2]   # x positions at last time step (even indices)
y_end = data[-1, 1::2]  # y positions at last time step (odd indices)
# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))    
# Plot starting positions
ax1.scatter(x_start, y_start, color='steelblue', alpha=0.7, s=10, edgecolor='black', linewidth=0.5)
ax1.set_xlabel('$x$ ', fontsize=13, fontweight='bold')
ax1.set_ylabel('$y$ ', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_aspect('equal', adjustable='box')
# Plot ending positions
ax2.scatter(x_end, y_end, color='coral', alpha=0.7, s=10, edgecolor='black', linewidth=0.5)
ax2.set_xlabel('$x$ ', fontsize=13, fontweight='bold')
ax2.set_ylabel('$y$ ', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_aspect('equal', adjustable='box')
# Improve layout
ax1.set_xticks([0,50,100])  # Adjust range and step as needed
ax1.set_yticks([0,50,100])  # Adjust range and
ax2.set_xticks([0,50,100])  # Adjust range and step as needed
ax2.set_yticks([0,50,100])  # Adjust range and
ax1.set_xlim(0,100)
ax2.set_xlim(0,100)
ax1.set_ylim(0,100)
ax2.set_ylim(0,100)
plt.tight_layout()
plt.savefig('brownian_positions.png', dpi=300, bbox_inches='tight')
plt.show()
