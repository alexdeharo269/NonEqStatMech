import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Set professional seaborn style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

# Read data from text file
data = np.loadtxt('P1_3-data.txt')
x = data[:, 0]
y = data[:, 1]

# Find where there are breaks (large jumps indicate disconnected trajectories)
# Calculate distances between consecutive points
dx = np.diff(x)
dy = np.diff(y)
distances = np.sqrt(dx**2 + dy**2)

# Identify breaks (you can adjust threshold based on your data)
threshold = np.median(distances) * 5  # 5x median distance
breaks = np.where(distances > threshold)[0] + 1
breaks = np.concatenate([[0], breaks, [len(x)]])

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Plot each continuous segment separately
colors = sns.color_palette("husl", len(breaks)-1)
for i in range(len(breaks)-1):
    start = breaks[i]
    end = breaks[i+1]
    
    ax.plot(x[start:end], y[start:end], 
            color=colors[i % len(colors)], 
            linewidth=1.5, 
            alpha=0.7,
            marker='o',
            markersize=2,
            markevery=max(1, (end-start)//20))  # Show some markers

# Customize the plot
ax.set_xlabel('$x$ ', fontsize=13, fontweight='bold')
ax.set_ylabel('$y$ ', fontsize=13, fontweight='bold')
ax.grid(True, linestyle='--', alpha=0.3)
ax.set_aspect('equal', adjustable='box')

# Add subtle background
ax.set_xticks([-1,0,1])  # Adjust range and step as needed
ax.set_yticks([-1,0,1])  # Adjust range and step as needed
ax.set_facecolor('#fafafa')

plt.tight_layout()
plt.savefig('brownian_trajectories.png', dpi=300, bbox_inches='tight')
plt.show()

