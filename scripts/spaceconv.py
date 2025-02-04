import numpy as np
import matplotlib.pyplot as plt
import os

# Example data: replace these with your actual calculated errors and grid sizes
# errors =[]
# #size =[]
# dir = os.fsencode("../resources/conv")
# for file in sorted(os.listdir(dir)):   
#     data = np.genfromtxt(os.fsdecode(os.path.join(dir,file)), delimiter=",")
#     #size.append(int(os.fsdecode(file).split(".")[0]))
#     errors.append(data[-1, 2])
    
# print(errors)
#print(size)
errors = [7.50053e-05, 1.20893e-5, 2.98394e-6]
grid_sizes = np.array([20 , 50 , 100])
errors = np.array(errors)  # Replace with actual error values

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
p1 = np.log(errors[1] / errors[0]) / np.log(grid_sizes[0] / grid_sizes[1])
p2 = np.log(errors[2] / errors[1]) / np.log(grid_sizes[1] / grid_sizes[2])
print(f'Estimated order of convergence 1: {p1:.2f}')
print(f'Estimated order of convergence 2: {p2:.2f}')