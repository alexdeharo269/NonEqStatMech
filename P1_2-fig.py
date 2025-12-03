import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set seaborn style for professional appearance
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.4)
sns.set_palette("husl")

# Read data from file
data = np.loadtxt('Gaussian_stats.txt')
x_values = data[:, 0]
y_values = data[:, 1]

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

# Plot histogram for x values with KDE
sns.histplot(x_values, bins=50, kde=True, ax=ax1, 
             color='steelblue', alpha=0.6, edgecolor='black', linewidth=0.5)
ax1.set_xlabel('$X$', fontsize=14, fontweight='bold')
ax1.set_ylabel('$P(X)$', fontsize=14, fontweight='bold')
ax1.set_title('$X$ Distribution', fontsize=16, fontweight='bold', pad=15)
ax1.grid(True, alpha=0.3, linestyle='--')

# Add statistics text box for X
mean_x = np.mean(x_values)
std_x = np.std(x_values)
textstr_x = f'$\mu$ = {mean_x:.3f}\n$\sigma$ = {std_x:.3f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
ax1.text(0.05, 0.95, textstr_x, transform=ax1.transAxes, fontsize=11,
         verticalalignment='top', bbox=props)

# Plot histogram for y values with KDE
sns.histplot(y_values, bins=50, kde=True, ax=ax2,
             color='coral', alpha=0.6, edgecolor='black', linewidth=0.5)
ax2.set_xlabel('$Y$', fontsize=14, fontweight='bold')
ax2.set_ylabel('$P(Y)$', fontsize=14, fontweight='bold')
ax2.set_title('$Y$ Distribution', fontsize=16, fontweight='bold', pad=15)
ax2.grid(True, alpha=0.3, linestyle='--')

# Add statistics text box for Y
mean_y = np.mean(y_values)
std_y = np.std(y_values)
textstr_y = f'$\mu$ = {mean_y:.3f}\n$\sigma$ = {std_y:.3f}'
ax2.text(0.05, 0.95, textstr_y, transform=ax2.transAxes, fontsize=11,
         verticalalignment='top', bbox=props)

# Improve layout
plt.tight_layout()
plt.savefig('gaussian_distributions.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("=" * 50)
print("STATISTICAL SUMMARY")
print("=" * 50)
print(f"X Distribution:")
print(f"  Mean (μ):     {mean_x:.6f}")
print(f"  Std Dev (σ):  {std_x:.6f}")
print(f"  Min:          {np.min(x_values):.6f}")
print(f"  Max:          {np.max(x_values):.6f}")
print(f"\nY Distribution:")
print(f"  Mean (μ):     {mean_y:.6f}")
print(f"  Std Dev (σ):  {std_y:.6f}")
print(f"  Min:          {np.min(y_values):.6f}")
print(f"  Max:          {np.max(y_values):.6f}")
print("=" * 50)