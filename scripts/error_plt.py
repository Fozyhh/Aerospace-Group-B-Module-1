import numpy as np
import matplotlib.pyplot as plt

# Example data: replace these with your actual calculated errors and grid sizes
grid_sizes = np.array([50, 100])
errors = np.array([0.000111959,4.66244e-05])  # Replace with actual error values

# Plot the convergence plot with respect to grid sizes
plt.figure(figsize=(8, 6))
plt.loglog(grid_sizes, errors, marker='o', label='Numerical Error', linewidth=2)
plt.xlabel('log(Number of Cells)')
plt.ylabel('log(Error)')
plt.title('Convergence Plot')

# Plot a second-order reference line
error_ref = errors[0]
reference_line = error_ref * (grid_sizes / grid_sizes[0])**(-2)  # Negative slope for 2nd order
plt.loglog(grid_sizes, reference_line, '--', label='Second-Order Reference', linewidth=1.5)

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Save the plot as a file
plt.savefig('convergence_plot.png', dpi=300)  # Save as a high-resolution PNG file

plt.show()

# Calculate the convergence rate between the first and last points
p = np.log(errors[1] / errors[0]) / np.log(grid_sizes[1] / grid_sizes[0])
print(f'Estimated order of convergence: {p:.2f}')
