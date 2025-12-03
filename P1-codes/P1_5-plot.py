import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set seaborn style for professional appearance
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.4)

# Read data from file with proper column names
data = pd.read_csv('P1_5-data.txt', delim_whitespace=True, comment='#', 
                   names=['Gamma', 'D_noPBC', 'D_PBC'])

# Create log-log plot
plt.figure(figsize=(8, 5))

# Plot both datasets
sns.scatterplot(data=data, x='Gamma', y='D_noPBC', 
                color='steelblue', alpha=0.9, s=80, 
                marker='x', linewidth=1.5,
                label='D (no PBC)')

sns.scatterplot(data=data, x='Gamma', y='D_PBC', 
                color='coral', alpha=0.9, s=80, 
                marker='+', linewidth=1.5,
                label='D (with PBC)')

plt.xscale('log')
plt.yscale('log')

# Calculate log-log linear fits for both
log_x = np.log10(data['Gamma'])

# Fit for D_noPBC
log_y_noPBC = np.log10(data['D_noPBC'])
slope_noPBC, intercept_noPBC = np.polyfit(log_x, log_y_noPBC, 1)

# Fit for D_PBC
log_y_PBC = np.log10(data['D_PBC'])
slope_PBC, intercept_PBC = np.polyfit(log_x, log_y_PBC, 1)

# Add fit lines
x_fit = np.logspace(np.min(log_x), np.max(log_x), 100)

y_fit_noPBC = 10**intercept_noPBC * x_fit**slope_noPBC
plt.plot(x_fit, y_fit_noPBC, 'b--', linewidth=1, alpha=0.7,
         label=f'No PBC fit: slope = {slope_noPBC:.3f}')

y_fit_PBC = 10**intercept_PBC * x_fit**slope_PBC
plt.plot(x_fit, y_fit_PBC, 'r--', linewidth=1, alpha=0.7, 
         label=f'PBC fit: slope = {slope_PBC:.3f}')


plt.xlabel('$\Gamma$', fontsize=14, fontweight='bold')
plt.ylabel('Diffusivity $D$', fontsize=14)
plt.grid(True, alpha=0.3, linestyle='--')

# Improve legend
plt.legend(fontsize=10, framealpha=0.95, loc='best')

# Improve layout and save
plt.tight_layout()
plt.savefig('loglog_plot_comparison.png', dpi=300, bbox_inches='tight')

print(f"\nFit results:")
print(f"No PBC: slope = {slope_noPBC:.4f}, intercept = {intercept_noPBC:.4f}")
print(f"With PBC: slope = {slope_PBC:.4f}, intercept = {intercept_PBC:.4f}")
print(f"\nBoth slopes should be close to 1.0 (Einstein relation: D = Γ/γ with γ=1)")

plt.show()