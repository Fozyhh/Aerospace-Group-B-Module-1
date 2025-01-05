import numpy as np
import matplotlib.pyplot as plt

# Example data: replace these with your actual calculated errors and grid sizes
grid_sizes = np.array([10, 20, 40])
errors = np.array([0.215167, 0.0606426, 0.0165213])  # Replace with actual error values

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
#p = np.log(errors[1] / errors[0]) / np.log(grid_sizes[1] / grid_sizes[0])
#print(f'Estimated order of convergence: {p:.2f}')

#compute the convergence rate for all steps and print the average
#p = np.log(errors[1:] / errors[:-1]) / np.log(grid_sizes[1:] / grid_sizes[:-1])
p = np.log(errors[:-1] / errors[1:]) / np.log(grid_sizes[1:] / grid_sizes[:-1])
print(f'Estimated order of convergence: {p}')
print(f'Average order of convergence: {np.mean(p):.2f}')
print(f'Average order of convergence: {np.mean(p):.2f} +/- {np.std(p):.2f}')
